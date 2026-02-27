'''
parsl_config.py

Parsl configuration for the PLATLAS locus-defining pipeline on Aurora.

Uses LocalProvider so that parsl distributes work across nodes already
allocated to the outer PBS job — no secondary PBS job submissions.  The outer
PBS script (run_locus_defining.sh) must therefore request all worker nodes
upfront (select=32).

Workers per node is set conservatively because the bottleneck is I/O:
reading and decompressing large bgzip summary statistics files (GWAMA output
can be several GB uncompressed; MVP sumstats are smaller but still I/O-bound).

Parallelism summary:
  workers_per_node (4) x num_nodes (read from PBS_NODEFILE at runtime) = concurrent jobs
  Default: 4 × 32 = 128 concurrent locus-defining runs

***Run with ALCF frameworks module: module load frameworks***
'''

import os
from parsl.config import Config
from parsl.providers import LocalProvider
from parsl.executors import HighThroughputExecutor
from parsl.launchers import MpiExecLauncher
from parsl.addresses import address_by_hostname

# Directory where parsl writes worker logs and working files
execute_dir = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs"

# Concurrent locus-defining jobs per node (conservative: large compressed files, heavy IO)
workers_per_node = 4

# Determine number of allocated nodes from PBS environment at import time.
# PBS_NODEFILE contains one hostname per line (one entry per node when mpiprocs=1).
_nodefile = os.environ.get("PBS_NODEFILE", "")
if _nodefile and os.path.exists(_nodefile):
    with open(_nodefile) as _fh:
        num_nodes = len(set(line.strip() for line in _fh if line.strip()))
else:
    num_nodes = 1  # fallback for interactive/testing use

locus_defining_config = Config(
    executors=[
        HighThroughputExecutor(
            # Address of the interchange on the head node; remote workers use
            # this to connect back.  Required when workers run on other nodes.
            address=address_by_hostname(),
            # Number of tasks that can run simultaneously on each node
            max_workers_per_node=workers_per_node,
            # Buffer extra tasks so workers stay busy
            prefetch_capacity=workers_per_node * 2,
            provider=LocalProvider(
                # Use all nodes already allocated to this PBS job
                nodes_per_block=num_nodes,
                # One parsl manager process per node; it dispatches to workers
                launcher=MpiExecLauncher(
                    bind_cmd="--cpu-bind", overrides="--ppn 1"
                ),
                # Environment setup run before each worker process launches.
                # NUMEXPR_MAX_THREADS is set to 208/4=52 so that all 4 workers
                # on a node together use at most 208 threads (matching the node
                # CPU count) and numexpr does not emit an oversubscription warning.
                worker_init=(
                    "source $HOME/.bashrc; "
                    "module load frameworks; "
                    "export NUMEXPR_MAX_THREADS=52; "
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
