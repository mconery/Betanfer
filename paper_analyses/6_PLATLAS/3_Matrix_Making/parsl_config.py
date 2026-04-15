'''
parsl_config.py

Parsl configuration for the PLATLAS LD matrix-making pipeline on Aurora.

Uses LocalProvider so that parsl distributes work across nodes already
allocated to the outer PBS job — no secondary PBS job submissions.  The outer
PBS script (run_matrix_making.sh) must therefore request all worker nodes
upfront.

Each Aurora GPU node has 6 GPUs × 2 tiles = 12 GPU tiles.  One parsl worker
is bound to one tile per physical GPU (6 workers/node) via:
  - available_accelerators: sets ZE_AFFINITY_MASK per worker (Intel XPU)
  - cpu_affinity:           pins each worker to its own 32 logical CPU cores
                            (16 physical + 16 HT siblings), preventing CPU
                            oversubscription across workers on the same node.

Using 6 workers instead of 12 halves the concurrent Lustre file operations per
node, preventing client-side RPC saturation when processing large cohorts.

betanfer infer is GPU-accelerated and detects its assigned tile automatically
from the ZE_AFFINITY_MASK environment variable set by parsl.  The --num-processes
option in the betanfer command controls parallel CPU use for the p-value
filtering and beta extraction stages; 12 processes uses 12 of the 16 physical
cores available in each worker's affinity group.

Parallelism summary:
  6 workers/node × num_nodes (from PBS_NODEFILE) = concurrent betanfer runs

***Run with the parsl_gpu conda environment: conda activate /flare/GeomicVar/mconery/tools/parsl_gpu***
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

# 6 workers per Aurora node: one per physical GPU (skipping the second tile of
# each GPU).  Aurora has 6 GPUs × 2 tiles = 12 tiles total; in FLAT hierarchy
# tiles appear as independent Level Zero devices 0–11.  Using every other tile
# (0, 2, 4, 6, 8, 10) halves the number of concurrent workers per node, which
# reduces Lustre client contention when processing large cohorts like FINNGEN
# that require thousands of concurrent file reads.
tile_names = [str(i) for i in range(0, 12, 2)]

# CPU affinity: one group per worker (16 physical + 16 HT siblings = 32 logical
# CPUs).  Adjacent 8-core groups from the original 12-worker layout are merged
# in pairs, maintaining NUMA locality (workers 0–2 → socket 0, 3–5 → socket 1).
# Total CPUs pinned: 6 × 32 = 192 (out of 208; remaining 16 reserved for OS).
cpu_affinity = (
    "list:"
    "1-16,105-120:"
    "17-32,121-136:"
    "33-48,137-152:"
    "53-68,157-172:"
    "69-84,173-188:"
    "85-100,189-204"
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
                # parsl_gpu conda env provides parsl, PyTorch with Intel XPU (IPEX)
                # support, pysam, and betanfer.
                worker_init=(
                    "source /flare/GeomicVar/mconery/tools/miniconda3/etc/profile.d/conda.sh; "
                    "conda activate /flare/GeomicVar/mconery/tools/parsl_gpu; "
                    "export ZE_FLAT_DEVICE_HIERARCHY=FLAT; "
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
