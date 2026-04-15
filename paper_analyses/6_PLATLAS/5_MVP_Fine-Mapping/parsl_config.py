'''
parsl_config.py

Parsl configuration for the PLATLAS MVP fine-mapping pipeline on Aurora.

Uses LocalProvider so that parsl distributes work across nodes already
allocated to the outer PBS job — no secondary PBS job submissions.  The outer
PBS script (run_mvp_fine_mapping.sh) must therefore request all worker nodes
upfront.

Workers per node is set to 4 (I/O and compute bound: each job tabix-extracts
a locus, runs prep_z-score_files.py, then runs a fine-mapping tool that reads
a large LD matrix).

Parallelism summary:
  workers_per_node (4) x num_nodes (read from PBS_NODEFILE at runtime) = concurrent jobs
  Default: 4 × 16 = 64 concurrent fine-mapping runs

Both the outer PBS script (run_mvp_fine_mapping.sh) and the worker_init below
activate the `finemap` conda environment.  This is required: parsl's interchange
(head node) and managers (worker nodes) must run the same Python version or they
cannot communicate.  The finemap env also provides tabix, RSparsePro (Python),
and SuSiE (Rscript) needed by the per-locus wrapper scripts.

***Run in Conda finemap environment***
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
# Aurora nodes have 104 CPUs each.  Setting MKL_NUM_THREADS / OMP_NUM_THREADS
# below to 104 // workers_per_node keeps BLAS from oversubscribing across workers.
workers_per_node = 8

# Determine number of allocated nodes from PBS environment at import time.
# PBS_NODEFILE contains one hostname per line (one entry per node when mpiprocs=1).
_nodefile = os.environ.get("PBS_NODEFILE", "")
if _nodefile and os.path.exists(_nodefile):
    with open(_nodefile) as _fh:
        num_nodes = len(set(line.strip() for line in _fh if line.strip()))
else:
    num_nodes = 1  # fallback for interactive/testing use

mvp_finemap_config = Config(
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
                # Activate the finemap conda env for RSparsePro (Python) and
                # SuSiE (Rscript).  Thread limits prevent numpy/MKL and R/OpenBLAS
                # from oversubscribing CPUs across co-resident workers.
                # With 8 workers on 104-CPU nodes: 104 // 8 = 13 threads each.
                worker_init=(
                    "source $HOME/.bashrc; "
                    "conda activate /flare/GeomicVar/mconery/tools/miniconda3/envs/finemap; "
                    "export NUMEXPR_MAX_THREADS=13; "
                    "export MKL_NUM_THREADS=13; "
                    "export OMP_NUM_THREADS=13; "
                    "export OPENBLAS_NUM_THREADS=13; "
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
