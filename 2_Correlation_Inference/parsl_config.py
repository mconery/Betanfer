'''
parsl_config.py

The purpose of this code is to hold parsl configurations needed for running 
different tasks on Aurora.

***Run in Conda finemap_old environment***
'''

# config.py
import os
from parsl.config import Config

# PBSPro is the right provider for ALCF:
from parsl.providers import PBSProProvider
# The high throughput executor is for scaling large single core/tile/gpu tasks on HPC system:
from parsl.executors import HighThroughputExecutor
# Use the MPI launcher to launch worker processes:
from parsl.launchers import MpiExecLauncher

# These options will run work in 1 node batch jobs run one at a time
nodes_per_job = 1
max_num_jobs = 1
tile_names = [f'{gid}.{tid}' for gid in range(6) for tid in range(2)]

# The config will launch workers from this directory
execute_dir = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs"

aurora_single_tile_config = Config(
    executors=[
        HighThroughputExecutor(
            # Ensures one worker per GPU tile on each node
            available_accelerators=tile_names,
            max_workers_per_node=208,
            # Distributes threads to workers/tiles in a way optimized for Aurora
            cpu_affinity="list:1-8,105-112:9-16,113-120:17-24,121-128:25-32,129-136:33-40,137-144:41-48,145-152:53-60,157-164:61-68,165-172:69-76,173-180:77-84,181-188:85-92,189-196:93-100,197-204",
            # Increase if you have many more tasks than workers
            prefetch_capacity=0,
            # Options that specify properties of PBS Jobs
            provider=PBSProProvider(
                # Project name
                account="GeomicVar",
                # Submission queue
                queue="debug",
                # Commands run before workers launched
                # Make sure to activate your environment where Parsl is installed
                worker_init=f'''source $HOME/.bashrc; conda activate finemap_old; cd {execute_dir}''',
                # Wall time for batch jobs
                walltime="1:00:00",
                # Change if data/modules located on other filesystem
                scheduler_options="#PBS -l filesystems=home:flare",
                # Ensures 1 manger per node; the manager will distribute work to its 12 workers, one per tile
                launcher=MpiExecLauncher(bind_cmd="--cpu-bind", overrides="--ppn 1"),
                # options added to #PBS -l select aside from ncpus
                select_options="",
                # Number of nodes per PBS job
                nodes_per_block=nodes_per_job,
                # Minimum number of concurrent PBS jobs running workflow
                min_blocks=0,
                # Maximum number of concurrent PBS jobs running workflow
                max_blocks=max_num_jobs,
                # Hardware threads per node
                cpus_per_node=208,
            ),
        ),
    ],
    # How many times to retry failed tasks
    # this is necessary if you have tasks that are interrupted by a PBS job ending
    # so that they will restart in the next job
    retries=1,
)