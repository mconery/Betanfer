#!/usr/bin/env python3

"""
Queue monitoring script for distributed SuSiE fine-mapping.

This script monitors progress of the work queue in real-time and provides
status updates on completion percentage, in-progress batches, and remaining work.

Usage:
    python monitor_queue.py --queue-dir /lustre/work_queue

This runs continuously until all work is complete.
"""

import json
import time
import argparse
import sys
from pathlib import Path
from typing import Dict
import curses

def get_queue_status(queue_dir: str) -> Dict:
    """
    Get current status of the work queue.
    
    Args:
        queue_dir: Root directory of work queue
    
    Returns:
        Dictionary with status information
    """
    queue_path = Path(queue_dir)
    
    try:
        work_count = len(list((queue_path / "work").glob("*.json")))
        in_progress_count = len(list((queue_path / "in_progress").glob("*.json")))
        completed_count = len(list((queue_path / "completed").glob("*.json")))
        
        state_file = queue_path / "queue_state.json"
        if state_file.exists():
            with open(state_file) as f:
                state = json.load(f)
        else:
            state = {'total_commands': 0}
        
        total_commands = state.get('total_commands', 0)
        total_batches = state.get('total_batches', 0)
        
        # Calculate derived statistics
        total_work_items = work_count + in_progress_count + completed_count
        
        if total_commands > 0:
            percent_done = 100.0 * completed_count / total_batches
            percent_commands = 100.0 * completed_count / total_batches  # Assumes roughly equal batch sizes
        else:
            percent_done = 0
            percent_commands = 0
        
        return {
            'work_remaining': work_count,
            'in_progress': in_progress_count,
            'completed': completed_count,
            'total_batches': total_batches,
            'total_commands': total_commands,
            'percent_done': percent_done,
            'is_complete': work_count == 0 and in_progress_count == 0,
            'timestamp': time.time(),
        }
    
    except Exception as e:
        return {
            'error': str(e),
            'work_remaining': -1,
            'in_progress': -1,
            'completed': -1,
            'total_batches': -1,
            'total_commands': -1,
            'percent_done': 0,
            'is_complete': False,
            'timestamp': time.time(),
        }


def format_status(status: Dict, elapsed: float) -> str:
    """
    Format status dictionary into human-readable string.
    
    Args:
        status: Status dictionary from get_queue_status
        elapsed: Seconds elapsed since monitoring started
    
    Returns:
        Formatted status string
    """
    if 'error' in status:
        return f"Error: {status['error']}"
    
    lines = []
    lines.append(f"Time elapsed: {int(elapsed):,d}s ({elapsed/60:.1f}m)")
    lines.append("")
    lines.append("BATCH STATUS:")
    lines.append(f"  Completed: {status['completed']:6d} / {status['total_batches']:6d} "
                f"({status['percent_done']:5.1f}%)")
    lines.append(f"  In Progress: {status['in_progress']:6d}")
    lines.append(f"  Remaining: {status['work_remaining']:6d}")
    lines.append("")
    
    if status['total_batches'] > 0 and status['completed'] > 0:
        # Calculate ETA
        batches_per_second = status['completed'] / elapsed
        remaining_batches = status['work_remaining']
        estimated_seconds_remaining = remaining_batches / batches_per_second if batches_per_second > 0 else 0
        
        lines.append("PERFORMANCE:")
        lines.append(f"  Throughput: {batches_per_second:.2f} batches/sec")
        if estimated_seconds_remaining > 0:
            lines.append(f"  Estimated time remaining: {int(estimated_seconds_remaining):,d}s "
                        f"({estimated_seconds_remaining/60:.1f}m)")
        lines.append("")
    
    if status['is_complete']:
        lines.append("✓ ALL WORK COMPLETE!")
    
    return "\n".join(lines)


def monitor_queue_simple(queue_dir: str, interval: int = 5) -> None:
    """
    Simple text-based queue monitoring (no curses).
    
    Args:
        queue_dir: Root directory of work queue
        interval: Update interval in seconds
    """
    queue_path = Path(queue_dir)
    
    if not queue_path.exists():
        print(f"Error: Queue directory not found: {queue_dir}")
        sys.exit(1)
    
    print(f"Monitoring queue at {queue_dir}")
    print(f"Update interval: {interval} seconds")
    print("=" * 70)
    print("")
    
    start_time = time.time()
    last_status = None
    
    try:
        while True:
            status = get_queue_status(queue_dir)
            elapsed = time.time() - start_time
            
            # Print header with timestamp
            print(f"\n[{time.strftime('%Y-%m-%d %H:%M:%S')}] Queue Status")
            print("-" * 70)
            
            # Print status
            print(format_status(status, elapsed))
            
            # Check if complete
            if status['is_complete']:
                print("=" * 70)
                print(f"Monitoring complete at {time.strftime('%Y-%m-%d %H:%M:%S')}")
                print(f"Total runtime: {int(elapsed):,d}s ({elapsed/60:.1f}m)")
                break
            
            # Wait before next update
            time.sleep(interval)
    
    except KeyboardInterrupt:
        print("\n\nMonitoring stopped by user")
        elapsed = time.time() - start_time
        print(f"Runtime: {int(elapsed):,d}s ({elapsed/60:.1f}m)")
        sys.exit(0)


def monitor_queue_curses(stdscr, queue_dir: str, interval: int = 5) -> None:
    """
    Curses-based queue monitoring with live updates.
    
    Args:
        stdscr: Curses window object
        queue_dir: Root directory of work queue
        interval: Update interval in seconds
    """
    queue_path = Path(queue_dir)
    
    if not queue_path.exists():
        print(f"Error: Queue directory not found: {queue_dir}")
        return
    
    # Configure curses
    curses.curs_set(0)  # Hide cursor
    stdscr.nodelay(True)  # Non-blocking input
    
    start_time = time.time()
    last_update = 0
    
    try:
        while True:
            # Non-blocking check for 'q' key to quit
            try:
                key = stdscr.getkey()
                if key.lower() == 'q':
                    return
            except:
                pass
            
            # Update every interval seconds
            now = time.time()
            if now - last_update < interval:
                time.sleep(0.1)
                continue
            
            last_update = now
            status = get_queue_status(queue_dir)
            elapsed = now - start_time
            
            # Clear screen
            stdscr.clear()
            
            # Print to curses window
            try:
                stdscr.addstr(0, 0, "=" * 70)
                stdscr.addstr(1, 0, "SuSiE Fine-Mapping Queue Monitor (Press 'q' to quit)")
                stdscr.addstr(2, 0, "=" * 70)
                
                row = 4
                for line in format_status(status, elapsed).split("\n"):
                    stdscr.addstr(row, 0, line)
                    row += 1
                
                stdscr.refresh()
                
                # Exit if complete
                if status['is_complete']:
                    stdscr.addstr(row + 2, 0, "Press any key to exit...")
                    stdscr.refresh()
                    stdscr.getch()
                    return
            
            except curses.error:
                pass  # Window too small, just continue


def main():
    """Main entry point with argument parsing."""
    parser = argparse.ArgumentParser(
        description='Monitor SuSiE fine-mapping work queue progress',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic monitoring
  python monitor_queue.py --queue-dir /lustre/work_queue
  
  # With custom update interval (10 seconds)
  python monitor_queue.py --queue-dir /lustre/work_queue --interval 10
  
  # Disable curses (text-only output)
  python monitor_queue.py --queue-dir /lustre/work_queue --no-curses
        """
    )
    
    parser.add_argument(
        '--queue-dir',
        required=True,
        help='Root directory of work queue'
    )
    parser.add_argument(
        '--interval',
        type=int,
        default=5,
        help='Update interval in seconds (default: 5)'
    )
    parser.add_argument(
        '--no-curses',
        action='store_true',
        help='Disable curses, use simple text output'
    )
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.interval < 1:
        print(f"Error: interval must be >= 1, got {args.interval}")
        return 1
    
    queue_path = Path(args.queue_dir)
    if not queue_path.exists():
        print(f"Error: Queue directory not found: {args.queue_dir}")
        return 1
    
    # Run monitoring
    if args.no_curses or not sys.stdout.isatty():
        # Use simple text output if curses disabled or not a TTY
        monitor_queue_simple(args.queue_dir, args.interval)
    else:
        # Try to use curses for nicer output
        try:
            curses.wrapper(monitor_queue_curses, args.queue_dir, args.interval)
        except KeyboardInterrupt:
            pass
        except Exception as e:
            print(f"Curses error: {e}, falling back to simple monitoring")
            monitor_queue_simple(args.queue_dir, args.interval)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
