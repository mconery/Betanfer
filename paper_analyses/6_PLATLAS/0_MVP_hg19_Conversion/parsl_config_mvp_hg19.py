'''
parsl_config_mvp_hg19.py

The purpose of this code is to hold the parsl configuration needed for running
convert_mvp_to_hg19.py across all MVP AFR/EUR summary statistics files on
Aurora.  Exports `mvp_hg19_config`.

Uses LocalProvider so that parsl distributes work across nodes already
allocated to the outer PBS job — no secondary PBS job submissions.  The outer
PBS script (run_mvp_hg19_conversion.sh) must therefore request all worker
nodes upfront (select=32).

convert_mvp_to_hg19.py is single-threaded (pandas I/O and merge do not use
multiple CPUs).  The annotation file is 511 MB gzip (~3–5 GB in-memory) and
each sumstats file adds ~5–15 GB.  With 1 TB RAM per Aurora node, 4 workers
per node is well within memory limits (~25 GB × 4 = ~100 GB peak).  Workers
2–4 on the same node also benefit from the OS page cache for the annotation
file, which is read by worker 1 first.

Parallelism summary:
  workers_per_node (4) x num_nodes (read from PBS_NODEFILE at runtime) = concurrent jobs
  Default: 4 × 32 = 128 concurrent MVP hg19 conversion runs

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

# 4 workers per node: script is single-threaded; 1 TB RAM supports ~4 concurrent
# workers at ~25 GB peak each; page cache benefits workers 2-4 for annotation read
workers_per_node = 4

# Determine number of allocated nodes from PBS environment at import time.
# PBS_NODEFILE contains one hostname per line (one entry per node when mpiprocs=1).
_nodefile = os.environ.get("PBS_NODEFILE", "")
if _nodefile and os.path.exists(_nodefile):
    with open(_nodefile) as _fh:
        num_nodes = len(set(line.strip() for line in _fh if line.strip()))
else:
    num_nodes = 1  # fallback for interactive/testing use

mvp_hg19_config = Config(
    executors=[
        HighThroughputExecutor(
            # Address of the interchange on the head node; remote workers use
            # this to connect back.  Required when workers run on other nodes.
            address=address_by_hostname(),
            # 1 task at a time per node to limit peak memory use
            max_workers_per_node=workers_per_node,
            # Buffer 1 extra task so workers stay busy
            prefetch_capacity=workers_per_node,
            provider=LocalProvider(
                # Use all nodes already allocated to this PBS job
                nodes_per_block=num_nodes,
                # One parsl manager process per node; it dispatches to workers
                launcher=MpiExecLauncher(
                    bind_cmd="--cpu-bind", overrides="--ppn 1"
                ),
                # Environment setup run before each worker process launches.
                # NUMEXPR_MAX_THREADS=52 (208 CPUs / 4 workers per node).
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
