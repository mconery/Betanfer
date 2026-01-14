#!/usr/bin/env python3

"""
Node Manager for Distributed SuSiE Fine-Mapping

This script runs on each compute node and pulls work from a distributed queue
stored on Lustre. It manages local execution of SuSiE commands using multiple
threads and handles atomic work claiming from the queue.

Usage:
    python node_manager.py --queue-dir /path/to/queue --node-id 0 --workers 28

Expected queue structure:
    queue_dir/
        ├── work/                    (unclaimed batches)
        │   ├── batch_000000.json
        │   └── batch_000001.json
        ├── in_progress/             (currently processing)
        ├── completed/               (finished batches)
        └── queue_state.json         (metadata)
"""

import os
import json
import time
import subprocess
import argparse
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock
import logging
from typing import Dict, List, Tuple, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - [Node %(node_id)d] - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
    ]
)


class NodeWorkQueue:
    """
    Thread-safe work queue manager for a single compute node.
    
    Handles atomic claiming of batches from a shared Lustre work queue
    and manages local execution of commands.
    """
    
    def __init__(self, queue_dir: str, node_id: int, workers: int = 28):
        """
        Initialize node work queue manager.
        
        Args:
            queue_dir: Path to root work queue directory
            node_id: Unique identifier for this node
            workers: Number of parallel worker threads
        """
        self.queue_dir = Path(queue_dir)
        self.node_id = node_id
        self.workers = workers
        self.claimed_batches = set()
        self.claim_lock = Lock()
        self.logger = logging.getLogger(__name__)
        self.logger = self.logger.getChild(str(node_id))
        
        # Set up environment for R subprocesses
        self.env = os.environ.copy()
        # Force single-threaded BLAS to prevent oversubscription
        self.env['OPENBLAS_NUM_THREADS'] = '1'
        self.env['MKL_NUM_THREADS'] = '1'
        self.env['OMP_NUM_THREADS'] = '1'
        self.env['BLIS_NUM_THREADS'] = '1'
        self.env['VECLIB_MAXIMUM_THREADS'] = '1'
        
        # Verify queue directories exist
        self._verify_queue_structure()
        
        # Initialize statistics
        self.stats = {
            'batches_processed': 0,
            'commands_executed': 0,
            'commands_successful': 0,
            'commands_failed': 0,
            'total_runtime_seconds': 0.0,
        }
    
    def _verify_queue_structure(self) -> None:
        """Verify that queue directory structure exists."""
        required_dirs = ['work', 'in_progress', 'completed', 'logs']
        
        for dirname in required_dirs:
            dirpath = self.queue_dir / dirname
            if not dirpath.exists():
                self.logger.error(f"Queue directory missing: {dirpath}")
                raise FileNotFoundError(f"Required directory not found: {dirpath}")
        
        # Check for queue_state.json
        state_file = self.queue_dir / "queue_state.json"
        if not state_file.exists():
            self.logger.error(f"Queue state file missing: {state_file}")
            raise FileNotFoundError(f"Required file not found: {state_file}")
        
        self.logger.info("Queue structure verified")
    
    def claim_batch(self) -> Optional[Path]:
        """
        Atomically claim the next available batch from the work queue.
        
        Uses atomic file rename (mv) on Lustre to ensure no two nodes
        claim the same batch simultaneously.
        
        Returns:
            Path to claimed batch file in in_progress/, or None if no work
        """
        work_dir = self.queue_dir / "work"
        in_progress_dir = self.queue_dir / "in_progress"
        
        with self.claim_lock:
            # List all unclaimed batches
            try:
                available = sorted([
                    f for f in work_dir.iterdir() 
                    if f.suffix == '.json'
                ])
            except FileNotFoundError:
                self.logger.warning("Work directory not accessible")
                return None
            
            if not available:
                return None
            
            # Claim first available via atomic rename
            batch_file = available[0]
            batch_name = batch_file.stem
            
            try:
                in_progress_file = in_progress_dir / f"{batch_name}.json"
                # Atomic rename - Lustre guarantees this is atomic at the metadata level
                batch_file.rename(in_progress_file)
                self.claimed_batches.add(batch_name)
                
                self.logger.info(f"Claimed batch: {batch_name}")
                return in_progress_file
            
            except FileExistsError:
                self.logger.warning(f"Batch {batch_name} already in progress (race condition)")
                return None
            except Exception as e:
                self.logger.error(f"Error claiming batch: {e}")
                return None
    
    def mark_complete(self, batch_file: Path) -> None:
        """
        Mark a batch as completed and move to completed directory.
        
        Args:
            batch_file: Path to batch file in in_progress/
        """
        completed_dir = self.queue_dir / "completed"
        batch_name = batch_file.stem
        
        try:
            completed_file = completed_dir / f"{batch_name}.json"
            batch_file.rename(completed_file)
            self.claimed_batches.discard(batch_name)
            self.logger.info(f"Marked complete: {batch_name}")
        except Exception as e:
            self.logger.error(f"Error marking batch complete: {e}")
    
    def execute_batch(self, batch_file: Path) -> Tuple[int, int, int]:
        """
        Execute all commands in a batch using parallel workers.
        
        Args:
            batch_file: Path to batch JSON file
        
        Returns:
            Tuple of (total_commands, successful, failed)
        """
        # Read batch file
        try:
            with open(batch_file) as f:
                commands = json.load(f)
        except Exception as e:
            self.logger.error(f"Error reading batch file {batch_file}: {e}")
            return 0, 0, 0
        
        batch_name = batch_file.stem
        self.logger.info(f"Executing batch {batch_name} with {len(commands)} commands")
        
        # Create directory for batch logs if needed
        log_dir = self.queue_dir / "logs"
        batch_log_file = log_dir / f"{batch_name}.log"
        
        # Execute commands in parallel
        successful = 0
        failed = 0
        start_time = time.time()
        
        with ThreadPoolExecutor(max_workers=self.workers) as executor:
            futures = [
                executor.submit(self._run_command, cmd_data, idx, batch_name)
                for idx, cmd_data in enumerate(commands)
            ]
            
            # Process results as they complete
            for future in as_completed(futures):
                try:
                    cmd, return_code, stderr = future.result()
                    if return_code == 0:
                        successful += 1
                    else:
                        failed += 1
                except Exception as e:
                    self.logger.error(f"Exception executing command: {e}")
                    failed += 1
        
        elapsed = time.time() - start_time
        
        # Log batch results
        self.logger.info(
            f"Batch {batch_name} complete: "
            f"{successful}/{len(commands)} successful, "
            f"{failed} failed in {elapsed:.2f}s"
        )
        
        # Write detailed log
        try:
            with open(batch_log_file, 'w') as f:
                f.write(f"Batch: {batch_name}\n")
                f.write(f"Node: {self.node_id}\n")
                f.write(f"Executed: {len(commands)} commands\n")
                f.write(f"Successful: {successful}\n")
                f.write(f"Failed: {failed}\n")
                f.write(f"Duration: {elapsed:.2f}s\n")
        except Exception as e:
            self.logger.warning(f"Could not write batch log: {e}")
        
        # Update statistics
        self.stats['batches_processed'] += 1
        self.stats['commands_executed'] += len(commands)
        self.stats['commands_successful'] += successful
        self.stats['commands_failed'] += failed
        self.stats['total_runtime_seconds'] += elapsed
        
        return len(commands), successful, failed
    
    def _run_command(self, cmd_data: Dict, idx: int, batch_name: str) -> Tuple[str, int, str]:
        """
        Execute a single command in a subprocess.
        
        Args:
            cmd_data: Dictionary with 'command' and 'metadata' keys
            idx: Command index within batch
            batch_name: Name of batch (for logging)
        
        Returns:
            Tuple of (command, return_code, stderr)
        """
        command = cmd_data.get('command', '')
        
        if not command:
            self.logger.error(f"Batch {batch_name} command {idx}: empty command")
            return '', -1, 'Empty command'
        
        try:
            result = subprocess.run(
                command,
                shell=True,
                capture_output=True,
                text=True,
                timeout=None,  # No timeout - let commands run as long as needed
                env=self.env
            )
            
            if result.returncode != 0:
                self.logger.debug(
                    f"Batch {batch_name} cmd {idx} failed (code {result.returncode}): "
                    f"{command[:100]}..."
                )
            
            return command, result.returncode, result.stderr[:500]
        
        except subprocess.TimeoutExpired:
            error_msg = "Command timed out"
            self.logger.error(f"Batch {batch_name} cmd {idx}: {error_msg}")
            return command, -1, error_msg
        
        except Exception as e:
            error_msg = f"Exception during execution: {str(e)}"
            self.logger.error(f"Batch {batch_name} cmd {idx}: {error_msg}")
            return command, -1, error_msg
    
    def log_statistics(self) -> None:
        """Log final statistics for this node."""
        self.logger.info("=== Final Statistics ===")
        self.logger.info(f"Batches processed: {self.stats['batches_processed']}")
        self.logger.info(f"Commands executed: {self.stats['commands_executed']}")
        self.logger.info(f"Commands successful: {self.stats['commands_successful']}")
        self.logger.info(f"Commands failed: {self.stats['commands_failed']}")
        self.logger.info(f"Total runtime: {self.stats['total_runtime_seconds']:.2f}s")
        
        if self.stats['commands_executed'] > 0:
            success_rate = 100.0 * self.stats['commands_successful'] / self.stats['commands_executed']
            avg_cmd_time = self.stats['total_runtime_seconds'] / self.stats['commands_executed']
            self.logger.info(f"Success rate: {success_rate:.1f}%")
            self.logger.info(f"Avg time per command: {avg_cmd_time:.2f}s")
    
    def worker_loop(self, idle_timeout: int = 10) -> None:
        """
        Main worker loop for this node.
        
        Continuously polls for available work and executes batches until
        the queue is empty. Exits gracefully when no work is available
        for a specified timeout period.
        
        Args:
            idle_timeout: Number of consecutive empty polls before exiting (each ~0.5s)
        """
        self.logger.info(f"Starting work queue polling (will exit after {idle_timeout*0.5:.1f}s idle)")
        
        loop_start = time.time()
        idle_count = 0
        
        while idle_count < idle_timeout:
            batch_file = self.claim_batch()
            
            if batch_file is None:
                idle_count += 1
                
                # Log progress periodically
                if idle_count == 1:
                    self.logger.info(
                        f"No work available (idle count {idle_count}/{idle_timeout})"
                    )
                elif idle_count % 4 == 0:
                    self.logger.info(
                        f"Still waiting for work ({idle_count}/{idle_timeout})"
                    )
                
                time.sleep(0.5)
                continue
            
            idle_count = 0
            
            try:
                total, successful, failed = self.execute_batch(batch_file)
                self.mark_complete(batch_file)
            
            except Exception as e:
                self.logger.error(f"Error processing batch: {e}")
                # Move back to work queue on error for potential retry
                try:
                    work_file = self.queue_dir / "work" / batch_file.name
                    batch_file.rename(work_file)
                    self.logger.info("Moved batch back to work queue for retry")
                except Exception as e2:
                    self.logger.error(f"Could not recover batch: {e2}")
        
        total_time = time.time() - loop_start
        self.logger.info(
            f"Shutting down (no work available for {idle_timeout*0.5:.1f}s, "
            f"total runtime {total_time:.1f}s)"
        )
        self.log_statistics()


def main():
    """Main entry point with argument parsing."""
    parser = argparse.ArgumentParser(
        description='Node manager for distributed SuSiE fine-mapping',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python node_manager.py --queue-dir /lustre/queue --node-id 0 --workers 28
  
  # With custom idle timeout
  python node_manager.py --queue-dir /lustre/queue --node-id 5 --workers 14 --idle-timeout 20
        """
    )
    
    parser.add_argument(
        '--queue-dir',
        required=True,
        help='Path to root work queue directory on Lustre'
    )
    parser.add_argument(
        '--node-id',
        type=int,
        required=True,
        help='Unique node identifier (0 to N-1)'
    )
    parser.add_argument(
        '--workers',
        type=int,
        default=28,
        help='Number of parallel worker threads (default: 28)'
    )
    parser.add_argument(
        '--idle-timeout',
        type=int,
        default=10,
        help='Polls until idle before exiting (default: 10 = ~5 seconds)'
    )
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.workers < 1:
        print(f"Error: workers must be >= 1, got {args.workers}")
        return 1
    
    if args.idle_timeout < 1:
        print(f"Error: idle-timeout must be >= 1, got {args.idle_timeout}")
        return 1
    
    try:
        # Create queue manager
        queue = NodeWorkQueue(
            args.queue_dir,
            args.node_id,
            args.workers
        )
        
        # Run main loop
        queue.worker_loop(idle_timeout=args.idle_timeout)
        
        return 0
    
    except FileNotFoundError as e:
        logging.error(f"Queue directory error: {e}")
        return 1
    except KeyboardInterrupt:
        logging.info("Interrupted by user")
        return 0
    except Exception as e:
        logging.error(f"Fatal error: {e}")
        return 1


if __name__ == '__main__':
    sys.exit(main())
