#!/usr/bin/env python3
"""
Example script that loads a C++ Tomator1D input file and runs a simulation.

This demonstrates compatibility with existing input files.
"""

import sys
import os
import argparse
from pathlib import Path

# Add parent directory to path if running as script
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from tomator_dolfinx import (
    run_simulation,
    load_input_file,
)


def main():
    """Run simulation from a JSON input file."""
    
    parser = argparse.ArgumentParser(
        description="Run Tomator-dolfinx simulation from C++ input file"
    )
    parser.add_argument(
        "input_file",
        type=str,
        help="Path to JSON input file"
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default=None,
        help="Output directory (default: TomatorResults/<input_name>)"
    )
    parser.add_argument(
        "-t", "--tend",
        type=float,
        default=None,
        help="End time in seconds (overrides input file)"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print verbose output"
    )
    
    args = parser.parse_args()
    
    # Check input file exists
    input_path = Path(args.input_file)
    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}")
        sys.exit(1)
    
    print("=" * 60)
    print("Tomator-dolfinx - Python FEM Plasma Transport Solver")
    print("=" * 60)
    print(f"\nInput file: {input_path}")
    
    # Load input file
    params = load_input_file(str(input_path))
    
    if args.verbose:
        print(f"\nLoaded parameters:")
        for key, value in sorted(params.items()):
            print(f"  {key}: {value}")
    
    # Override end time if specified
    if args.tend is not None:
        params['tend'] = args.tend
    
    # Set output directory
    if args.output is not None:
        output_dir = args.output
    else:
        output_dir = f"TomatorResults/{input_path.stem}"
    
    print(f"Output directory: {output_dir}")
    print(f"Simulation time: {params.get('tend', 1e-3):.3e} s")
    
    # Run simulation
    print("\n" + "=" * 60)
    print("Starting simulation...")
    print("=" * 60 + "\n")
    
    state = run_simulation(params, output_dir=output_dir)
    
    print("\n" + "=" * 60)
    print("Simulation completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
