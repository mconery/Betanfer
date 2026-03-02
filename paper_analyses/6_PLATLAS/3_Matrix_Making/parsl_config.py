'''
parsl_config.py

Parsl configuration for the PLATLAS LD matrix-making pipeline on Aurora.

Uses LocalProvider so that parsl distributes work across nodes already
allocated to the outer PBS job — no secondary PBS job submissions.  The outer
PBS script (run_matrix_making.sh) must therefore request all worker nodes
upfront.

Each Aurora GPU node has 6 GPUs × 2 tiles = 12 GPU tiles.  One parsl worker
is bound to each tile via:
  - available_accelerators: sets ZE_AFFINITY_MASK per worker (Intel XPU)
  - cpu_affinity:           pins each worker to its own 16 logical CPU cores
                            (8 physical + 8 HT siblings), preventing CPU
                            oversubscription across workers on the same node.

betanfer infer is GPU-accelerated and detects its assigned tile automatically
from the ZE_AFFINITY_MASK environment variable set by parsl.  The --num-processes
option in the betanfer command controls parallel CPU use for the p-value
filtering stage; 8 processes matches the 8 physical cores per tile affinity
group.

Parallelism summary:
  12 workers/node × num_nodes (from PBS_NODEFILE) = concurrent betanfer runs

***Run with ALCF frameworks module: module load frameworks***
'''

import os
from parsl.config import Config
from parsl.providers import LocalProvider
from parsl.executors import HighThroughputExecutor
from parsl.launchers import MpiExecLauncher
from parsl.addresses import address_by_hostname

# Directory where parsl writes worker logs and working files
LOG_DIR = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs"

# Root of the betanfer repository (where betanfer is installed)
BETANFER_DIR = "/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer"

# 12 GPU tiles per Aurora node: 6 GPUs × 2 tiles each
# Tile names follow the ZE_AFFINITY_MASK convention "gpu_id.tile_id"
tile_names = [f'{gid}.{tid}' for gid in range(6) for tid in range(2)]

# CPU affinity: one group per tile (8 physical + 8 HT siblings = 16 logical CPUs).
# Source: https://docs.alcf.anl.gov/aurora/workflows/parsl/
# The 12 groups map 1:1 with tile_names above.  Total CPUs pinned: 12 × 16 = 192
# (out of 208 on the node; remaining 16 are reserved for OS / management).
cpu_affinity = (
    "list:"
    "1-8,105-112:"
    "9-16,113-120:"
    "17-24,121-128:"
    "25-32,129-136:"
    "33-40,137-144:"
    "41-48,145-152:"
    "53-60,157-164:"
    "61-68,165-172:"
    "69-76,173-180:"
    "77-84,181-188:"
    "85-92,189-196:"
    "93-100,197-204"
)

# Determine number of allocated nodes from the PBS environment.
# PBS_NODEFILE lists one hostname per allocated slot (one per node when mpiprocs=1).
_nodefile = os.environ.get("PBS_NODEFILE", "")
if _nodefile and os.path.exists(_nodefile):
    with open(_nodefile) as _fh:
        num_nodes = len(set(line.strip() for line in _fh if line.strip()))
else:
    num_nodes = 1  # fallback for interactive / testing use

matrix_making_config = Config(
    executors=[
        HighThroughputExecutor(
            # Address of the interchange on the head node; remote workers connect here
            address=address_by_hostname(),
            # Bind one worker per GPU tile
            available_accelerators=tile_names,
            # Pin each worker to its own CPU cores (prevents cross-worker contention)
            cpu_affinity=cpu_affinity,
            provider=LocalProvider(
                # All nodes already allocated to this PBS job
                nodes_per_block=num_nodes,
                # One parsl manager process per node dispatches to tile workers
                launcher=MpiExecLauncher(
                    bind_cmd="--cpu-bind", overrides="--ppn 1"
                ),
                # Environment setup run before each worker process launches.
                # frameworks module provides PyTorch with Intel XPU (IPEX) support
                # and pysam / betanfer must be installed before workers start.
                worker_init=(
                    "source $HOME/.bashrc; "
                    "module load frameworks; "
                    f"cd {BETANFER_DIR}"
                ),
                # Launch one block immediately; no dynamic scaling
                init_blocks=1,
                max_blocks=1,
            ),
        ),
    ],
    # Retry once if a task is interrupted (e.g. job resubmitted after walltime)
    retries=1,
)
