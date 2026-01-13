#!/usr/bin/env python3
"""
Flux Batch Executor: Execute commands from a batch file in parallel on a single node.

This script reads a command file (one command per line) and executes them using
a thread pool. It's designed to work with Flux's task allocation system, allowing
Flux to distribute batch files across nodes while individual nodes handle local
parallelization.

Usage:
    flux run -n 1 -c 32 python flux_batch_executor.py batch_commands_001.txt
    flux run -n 1 -c 32 python flux_batch_executor.py --workers 32 batch_commands_001.txt
"""

import argparse
import logging
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Tuple
import time


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("flux_batch_executor.log", mode="a"),
    ],
)
logger = logging.getLogger(__name__)


def execute_command(command: str, cmd_index: int) -> Tuple[int, str, int, str]:
    """
    Execute a single command in a subprocess.

    Args:
        command: Shell command to execute
        cmd_index: Index of this command in the batch

    Returns:
        Tuple of (cmd_index, command, return_code, stderr_output)
    """
    try:
        result = subprocess.run(
            command,
            shell=True,
            capture_output=True,
            text=True,
            timeout=None,  # No timeout; let commands run as long as needed
        )
        return (cmd_index, command, result.returncode, result.stderr)
    except subprocess.TimeoutExpired:
        error_msg = f"Command timed out"
        logger.error(f"[CMD {cmd_index}] {error_msg}: {command}")
        return (cmd_index, command, -1, error_msg)
    except Exception as e:
        error_msg = f"Exception during execution: {str(e)}"
        logger.error(f"[CMD {cmd_index}] {error_msg}: {command}")
        return (cmd_index, command, -1, error_msg)


def read_commands(batch_file: Path) -> List[str]:
    """
    Read commands from a batch file.

    Args:
        batch_file: Path to file containing commands (one per line)

    Returns:
        List of non-empty, non-comment commands
    """
    commands = []
    try:
        with open(batch_file, "r") as f:
            for line_num, line in enumerate(f, 1):
                # Strip whitespace
                line = line.strip()

                # Skip empty lines and comments
                if not line or line.startswith("#"):
                    continue

                commands.append(line)

        logger.info(f"Loaded {len(commands)} commands from {batch_file}")
        return commands

    except FileNotFoundError:
        logger.error(f"Batch file not found: {batch_file}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error reading batch file {batch_file}: {e}")
        sys.exit(1)


def execute_batch(
    commands: List[str],
    max_workers: int,
    report_interval: int = 50,
) -> Tuple[int, int, int]:
    """
    Execute all commands using a thread pool.

    Args:
        commands: List of commands to execute
        max_workers: Maximum number of parallel workers
        report_interval: Print progress every N completed commands

    Returns:
        Tuple of (total_commands, successful, failed)
    """
    total = len(commands)
    successful = 0
    failed = 0
    failed_commands = []

    logger.info(f"Starting execution of {total} commands with {max_workers} workers")
    start_time = time.time()

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks and create a mapping of futures to command indices
        futures = {
            executor.submit(execute_command, cmd, idx): idx
            for idx, cmd in enumerate(commands)
        }

        # Process results as they complete
        completed = 0
        for future in as_completed(futures):
            try:
                cmd_idx, cmd, return_code, stderr = future.result()
                completed += 1

                if return_code == 0:
                    successful += 1
                    if completed % report_interval == 0 or completed == total:
                        elapsed = time.time() - start_time
                        rate = completed / elapsed if elapsed > 0 else 0
                        logger.info(
                            f"Progress: {completed}/{total} completed "
                            f"({successful} OK, {failed} FAIL) - "
                            f"{rate:.1f} cmd/sec"
                        )
                else:
                    failed += 1
                    failed_commands.append((cmd_idx, cmd, return_code, stderr))
                    logger.warning(
                        f"[CMD {cmd_idx}] Failed with return code {return_code}\n"
                        f"  Command: {cmd}\n"
                        f"  Error: {stderr[:200]}"  # First 200 chars of error
                    )

            except Exception as e:
                failed += 1
                logger.error(f"Exception processing future: {e}")

    elapsed = time.time() - start_time
    logger.info(
        f"Batch execution completed in {elapsed:.2f}s - "
        f"Total: {total}, Successful: {successful}, Failed: {failed}"
    )

    # Report failures
    if failed_commands:
        logger.warning(f"\n{'='*60}")
        logger.warning(f"Summary of {failed} failed commands:")
        logger.warning(f"{'='*60}")
        for cmd_idx, cmd, return_code, stderr in failed_commands[:10]:  # Show first 10
            logger.warning(
                f"[CMD {cmd_idx}] Return code {return_code}: {cmd}\n"
                f"  Error: {stderr[:100]}"
            )
        if len(failed_commands) > 10:
            logger.warning(f"... and {len(failed_commands) - 10} more failures")

    return total, successful, failed


def main():
    parser = argparse.ArgumentParser(
        description="Execute commands from a batch file in parallel",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use default 32 workers
  python flux_batch_executor.py batch_commands.txt

  # Specify custom number of workers
  python flux_batch_executor.py --workers 56 batch_commands.txt

  # Set log level for debugging
  python flux_batch_executor.py --log-level DEBUG batch_commands.txt
        """,
    )

    parser.add_argument(
        "batch_file",
        type=Path,
        help="Path to batch file containing commands (one per line)",
    )
    parser.add_argument(
        "-w",
        "--workers",
        type=int,
        default=32,
        help="Number of parallel workers (default: 32)",
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging level (default: INFO)",
    )
    parser.add_argument(
        "--report-interval",
        type=int,
        default=50,
        help="Print progress every N commands (default: 50)",
    )

    args = parser.parse_args()

    # Adjust logging level
    logger.setLevel(getattr(logging, args.log_level))
    for handler in logger.handlers:
        handler.setLevel(getattr(logging, args.log_level))

    logger.info(f"Flux Batch Executor started")
    logger.info(f"  Batch file: {args.batch_file}")
    logger.info(f"  Max workers: {args.workers}")

    # Read commands
    commands = read_commands(args.batch_file)
    if not commands:
        logger.error("No commands found in batch file")
        sys.exit(1)

    # Execute batch
    total, successful, failed = execute_batch(
        commands,
        max_workers=args.workers,
        report_interval=args.report_interval,
    )

    # Exit with appropriate code
    if failed > 0:
        logger.error(f"Execution completed with {failed} failures")
        sys.exit(1)
    else:
        logger.info(f"All {total} commands completed successfully")
        sys.exit(0)


if __name__ == "__main__":
    main()
