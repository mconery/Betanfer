'''
parsl_config.py

Parsl configuration for the PLATLAS GWAS fine-mapping pipeline on Aurora.

Uses LocalProvider so that parsl distributes work across nodes already
allocated to the outer PBS job — no secondary PBS job submissions.  The outer
PBS script (run_gwas_fine_mapping.sh) must therefore request all worker nodes
upfront.

Workers per node is set to 4 (I/O and compute bound: each job tabix-extracts
per-cohort sumstats, runs a Python fine-mapping wrapper, and the fine-mapping
tool reads large per-cohort LD matrices).

Parallelism summary:
  workers_per_node (4) x num_nodes (read from PBS_NODEFILE at runtime) = concurrent jobs

Unlike the MVP fine-mapping step, this pipeline runs TWO different fine-mapping
methods (SuSiEx and MultiSuSiE) that require DIFFERENT conda environments.
The worker_init does NOT activate any specific conda env; instead, each command
in the command file uses the full path to its method's Python interpreter:
  SuSiEx:     /flare/GeomicVar/mconery/tools/miniconda3/envs/finemap/bin/python
  MultiSuSiE: /flare/GeomicVar/mconery/tools/miniconda3/envs/MultiSuSiE/bin/python

The outer PBS script runs under the finemap conda env (parsl is installed there).
Workers only need conda on the PATH, which `source $HOME/.bashrc` provides.

***Run with: conda activate finemap (parsl is installed there)***
'''

import os
from parsl.config import Config
from parsl.providers import LocalProvider
from parsl.executors import HighThroughputExecutor
from parsl.launchers import MpiExecLauncher
from parsl.addresses import address_by_hostname

# Directory where parsl writes worker logs and working files
execute_dir = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs"

# Concurrent fine-mapping jobs per node.
# Aurora nodes have 208 CPUs each.  Setting thread limits below to
# 208 // workers_per_node keeps numpy/MKL from oversubscribing across workers.
workers_per_node = 4

# Determine number of allocated nodes from PBS environment at import time.
# PBS_NODEFILE contains one hostname per line (one entry per node when mpiprocs=1).
_nodefile = os.environ.get("PBS_NODEFILE", "")
if _nodefile and os.path.exists(_nodefile):
    with open(_nodefile) as _fh:
        num_nodes = len(set(line.strip() for line in _fh if line.strip()))
else:
    num_nodes = 1  # fallback for interactive/testing use

gwas_finemap_config = Config(
    executors=[
        HighThroughputExecutor(
            # Address of the interchange on the head node; remote workers use
            # this to connect back.  Required when workers run on other nodes.
            address=address_by_hostname(),
            # Number of tasks that can run simultaneously on each node
            max_workers_per_node=workers_per_node,
            # Buffer extra tasks so workers don't go idle between completions
            prefetch_capacity=workers_per_node,
            provider=LocalProvider(
                # Use all nodes already allocated to this PBS job
                nodes_per_block=num_nodes,
                # One parsl manager process per node; it dispatches to workers
                launcher=MpiExecLauncher(
                    bind_cmd="--cpu-bind", overrides="--ppn 1"
                ),
                # Environment setup run before each worker process launches.
                # Activate finemap so parsl's process_worker_pool.py and Manager
                # run under the conda Python (system Python hits AF_UNIX path-too-long
                # when creating the multiprocessing Manager socket).  The actual
                # fine-mapping commands are unaffected: they embed full Python paths
                # for their own envs (finemap for SuSiEx, MultiSuSiE for MultiSuSiE).
                # Thread limits prevent numpy/MKL from oversubscribing CPUs across
                # co-resident workers.  With 4 workers on 208-CPU nodes: 52 threads.
                worker_init=(
                    "source $HOME/.bashrc; "
                    "conda activate /flare/GeomicVar/mconery/tools/miniconda3/envs/finemap; "
                    "export NUMEXPR_MAX_THREADS=52; "
                    "export MKL_NUM_THREADS=52; "
                    "export OMP_NUM_THREADS=52; "
                    "export OPENBLAS_NUM_THREADS=52; "
                    f"cd {execute_dir}"
                ),
                # Launch one block immediately; no dynamic scaling
                init_blocks=1,
                max_blocks=1,
            ),
        ),
    ],
    # Retry tasks once if interrupted (e.g. job resubmitted after walltime)
    retries=1,
)
