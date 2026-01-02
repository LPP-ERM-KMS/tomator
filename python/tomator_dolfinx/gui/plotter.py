#!/usr/bin/env python3
"""
Bokeh-based interactive plotter for Tomator simulation results.

This module provides functions to launch and stop the Bokeh visualization
server for real-time monitoring of simulation results.

Usage:
    # From Python code:
    from tomator_dolfinx.gui import launch_plotter, stop_plotter
    proc = launch_plotter("path/to/results.csv")
    # ... run simulation ...
    stop_plotter(proc)
    
    # From command line:
    python -m tomator_dolfinx.gui.plotter path/to/results.csv
"""

import os
import sys
import signal
import socket
import subprocess
import time
from pathlib import Path
from typing import Optional


# Global reference to plotter process for cleanup
_plotter_process: Optional[subprocess.Popen] = None


def find_free_port() -> int:
    """Find an available port for the Bokeh server."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("", 0))
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        return s.getsockname()[1]


def launch_plotter(csv_file: str, port: Optional[int] = None, show: bool = True) -> subprocess.Popen:
    """
    Launch the Bokeh plotter server for a CSV results file.
    
    Parameters
    ----------
    csv_file : str
        Path to the CSV results file to visualize.
    port : int, optional
        Port for the Bokeh server. If None, finds a free port automatically.
    show : bool
        If True, automatically opens the browser.
    
    Returns
    -------
    subprocess.Popen
        The Bokeh server process handle.
    """
    global _plotter_process
    
    csv_path = Path(csv_file).resolve()
    if not csv_path.parent.exists():
        csv_path.parent.mkdir(parents=True, exist_ok=True)
    
    if port is None:
        port = find_free_port()
    
    # Set environment variable for the plot_app to read
    env = os.environ.copy()
    env['TOMATOR_CSV_FILE'] = str(csv_path)
    
    # Path to plot_app.py
    plot_app_path = Path(__file__).parent / "plot_app.py"
    
    # Build bokeh serve command
    cmd = [
        sys.executable, "-m", "bokeh", "serve",
        str(plot_app_path),
        "--port", str(port),
    ]
    if show:
        cmd.append("--show")
    
    # Launch Bokeh server
    process = subprocess.Popen(
        cmd,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    
    _plotter_process = process
    
    print(f"[Plotter] Bokeh server started on port {port} (PID: {process.pid})")
    if show:
        print(f"[Plotter] Opening browser at http://localhost:{port}/plot_app")
    
    # Give Bokeh a moment to start
    time.sleep(1.0)
    
    return process


def stop_plotter(process: Optional[subprocess.Popen] = None) -> None:
    """
    Stop the Bokeh plotter server.
    
    Parameters
    ----------
    process : subprocess.Popen, optional
        The Bokeh server process to stop. If None, stops the global process.
    """
    global _plotter_process
    
    proc = process or _plotter_process
    
    if proc is None:
        return
    
    if proc.poll() is None:  # Process is still running
        print(f"[Plotter] Stopping Bokeh server (PID: {proc.pid})...")
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()
        print("[Plotter] Bokeh server stopped.")
    
    if proc is _plotter_process:
        _plotter_process = None


def _cleanup_on_exit():
    """Cleanup handler to stop plotter on program exit."""
    stop_plotter()


# Register cleanup handler
import atexit
atexit.register(_cleanup_on_exit)


def main():
    """Command-line entry point for the plotter."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Launch Bokeh plotter for Tomator simulation results"
    )
    parser.add_argument(
        "csv_file",
        type=str,
        help="Path to the CSV results file"
    )
    parser.add_argument(
        "-p", "--port",
        type=int,
        default=None,
        help="Port for Bokeh server (default: auto)"
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Don't automatically open browser"
    )
    
    args = parser.parse_args()
    
    csv_path = Path(args.csv_file)
    if not csv_path.exists():
        print(f"Error: CSV file not found: {csv_path}")
        print("Note: The plotter will wait for data to appear.")
    
    process = launch_plotter(
        str(csv_path),
        port=args.port,
        show=not args.no_show
    )
    
    print("\nPress Ctrl+C to stop the plotter...")
    
    try:
        # Wait for the process to finish or Ctrl+C
        process.wait()
    except KeyboardInterrupt:
        print("\nInterrupted by user.")
        stop_plotter(process)


if __name__ == "__main__":
    main()
