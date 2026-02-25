'''
parsl_config.py

The purpose of this code is to hold the parsl configuration needed for running
GWAMA fixed-effects meta-analyses on Aurora. Optimized for CPU-bound tasks
running multiple jobs per node (as opposed to GPU tile work).

Each PBS job requested by parsl runs on nodes_per_block nodes, with
workers_per_node concurrent GWAMA jobs per node. Workers are conservative
relative to CPU count because the bottleneck is I/O: decompressing large bgzip
summary statistics files (avg 1–1.8 GB compressed, est. 5–15 GB uncompressed)
before feeding them to GWAMA.

Parallelism summary:
  workers_per_node (4) x nodes_per_block (8) x max_blocks (1) = 32 concurrent jobs

The outer submission script (run_gwama_analysis.sh) requests 1 node to run the
parsl allocator; parsl then submits its own PBS worker jobs via PBSProProvider.

***Run in Conda finemap_old environment***
'''

import os
from parsl.config import Config
from parsl.providers import PBSProProvider
from parsl.executors import HighThroughputExecutor
from parsl.launchers import MpiExecLauncher

# Directory where parsl writes worker logs and working files
execute_dir = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs"

# Parallelism parameters — adjust if I/O or disk pressure warrants it
workers_per_node = 4   # concurrent GWAMA jobs per node
nodes_per_block  = 8   # nodes per PBS job submitted by parsl
max_num_blocks   = 1   # max concurrent PBS worker jobs (resubmit to continue)

gwama_meta_analysis_config = Config(
    executors=[
        HighThroughputExecutor(
            # Number of tasks that can run simultaneously on each node
            max_workers_per_node=workers_per_node,
            # Buffer extra tasks so workers stay busy between submissions
            prefetch_capacity=workers_per_node * 2,
            provider=PBSProProvider(
                # ALCF allocation account
                account="GeomicVar",
                # Aurora queue — capacity supports multi-node jobs
                queue="capacity",
                # Environment setup run before each worker process launches
                worker_init=(
                    f"source $HOME/.bashrc; "
                    f"conda activate finemap_old; "
                    f"cd {execute_dir}"
                ),
                # Wall time for each PBS worker job submitted by parsl
                walltime="48:00:00",
                # Lustre filesystem declarations required on Aurora
                scheduler_options="#PBS -l filesystems=home:flare",
                # One parsl manager process per node; it dispatches to workers
                launcher=MpiExecLauncher(
                    bind_cmd="--cpu-bind", overrides="--ppn 1"
                ),
                # No extra select options beyond ncpus (set by HighThroughputExecutor)
                select_options="",
                # Nodes per PBS job block
                nodes_per_block=nodes_per_block,
                # Scale down to 0 when idle
                min_blocks=0,
                # Maximum concurrent PBS worker jobs
                max_blocks=max_num_blocks,
                # Aurora node CPU count (used by parsl for resource accounting)
                cpus_per_node=208,
            ),
        ),
    ],
    # Retry tasks once if interrupted by a PBS job ending
    retries=1,
)
