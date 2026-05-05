'''
parsl_config.py

Parsl configuration for the 1000G LD matrix-making pipeline on Aurora.

Uses LocalProvider so that parsl distributes work across nodes already
allocated to the outer PBS job — no secondary PBS job submissions.  The outer
PBS script (run_1000g_matrix_making.sh) must therefore request all worker
nodes upfront.

Unlike the GPU-based betanfer matrix-making step, plink2 is natively
multi-threaded and uses all available CPUs within a single process.  One
parsl worker is therefore assigned per node, giving each plink2 job exclusive
use of all 208 logical CPUs on that node.

Parallelism summary:
  1 worker/node × num_nodes (from PBS_NODEFILE at runtime) = concurrent plink2 runs

***Run with: conda activate /flare/GeomicVar/mconery/tools/miniconda3/envs/finemap***
'''

import os
from parsl.config import Config
from parsl.providers import LocalProvider
from parsl.executors import HighThroughputExecutor
from parsl.launchers import MpiExecLauncher
from parsl.addresses import address_by_hostname

# Directory where parsl writes worker logs and working files
execute_dir = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs"

# One plink2 job per node; plink2 uses all 208 CPUs internally via --threads.
workers_per_node = 1

# Determine number of allocated nodes from PBS environment at import time.
# PBS_NODEFILE contains one hostname per line (one entry per node when mpiprocs=1).
_nodefile = os.environ.get("PBS_NODEFILE", "")
if _nodefile and os.path.exists(_nodefile):
    with open(_nodefile) as _fh:
        num_nodes = len(set(line.strip() for line in _fh if line.strip()))
else:
    num_nodes = 1  # fallback for interactive/testing use

matrix_1000g_config = Config(
    executors=[
        HighThroughputExecutor(
            # Address of the interchange on the head node; remote workers use
            # this to connect back.  Required when workers run on other nodes.
            address=address_by_hostname(),
            # One plink2 job per node (plink2 is multi-threaded)
            max_workers_per_node=workers_per_node,
            # Buffer one extra task so the worker doesn't idle between completions
            prefetch_capacity=workers_per_node,
            provider=LocalProvider(
                # Use all nodes already allocated to this PBS job
                nodes_per_block=num_nodes,
                # One parsl manager process per node; it dispatches to the worker
                launcher=MpiExecLauncher(
                    bind_cmd="--cpu-bind", overrides="--ppn 1"
                ),
                # Environment setup run before each worker process launches.
                # Activate finemap so parsl's process_worker_pool.py runs under
                # the conda Python.  Thread-limit variables prevent any incidental
                # numpy/MKL usage from over-subscribing CPUs; plink2 itself is
                # controlled via --threads in the command.
                worker_init=(
                    "source $HOME/.bashrc; "
                    "conda activate /flare/GeomicVar/mconery/tools/miniconda3/envs/finemap; "
                    "export NUMEXPR_MAX_THREADS=208; "
                    "export MKL_NUM_THREADS=208; "
                    "export OMP_NUM_THREADS=208; "
                    "export OPENBLAS_NUM_THREADS=208; "
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
