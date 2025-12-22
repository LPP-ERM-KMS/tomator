#!/usr/bin/env python3
"""
Debug script to diagnose blowup issues.
Runs just a few steps with detailed debug output.
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import numpy as np
from tomator_dolfinx.io.json_input import load_input_file
from tomator_dolfinx.mesh import create_mesh_from_geometry
from tomator_dolfinx.species import PlasmaState
from tomator_dolfinx.boundary import BoundaryConditions
from tomator_dolfinx.transport import TransportManager
from tomator_dolfinx.solver import BDF2Solver


def main():
    # Load TCV input file
    input_file = Path(__file__).parent / "TCV5151X_fixneDV_noH2.json"
    if not input_file.exists():
        input_file = Path(__file__).parent.parent.parent.parent / "examples/InputFiles/TCV5151X_fixneDV.json"
    
    print(f"Loading: {input_file}")
    params = load_input_file(str(input_file))
    
    # Print key params
    print("\nKey parameters:")
    print(f"  R = {params.get('R', '?')} m")
    print(f"  lHFS = {params.get('lHFS', '?')} m")
    print(f"  lLFS = {params.get('lLFS', '?')} m")
    print(f"  Dfix = {params.get('transport', {}).get('Dfix', '?')} m²/s")
    print(f"  Vfix = {params.get('transport', {}).get('Vfix', '?')} m/s")
    print(f"  dtinit = {params.get('time_step', {}).get('dtinit', '?')} s")
    print(f"  nevac = {params.get('initial_conditions', {}).get('nevac', '?')} m⁻³")
    
    # Create mesh
    mesh, radial_positions = create_mesh_from_geometry(
        R=params.get('R', 1.0),
        a=params.get('a', 0.2),
        lHFS=params.get('lHFS', 0.2),
        lLFS=params.get('lLFS', 0.2),
        num_cells=params.get('nmeshp', 100) - 1
    )
    
    print(f"\nMesh: r_min={radial_positions[0]:.4f} m, r_max={radial_positions[-1]:.4f} m, cells={len(radial_positions)-1}")
    
    # Create plasma state
    state = PlasmaState(mesh)
    
    # Add helium if enabled
    if params.get('bHe', False):
        state.add_helium_species()
    
    # Initialize
    state.initialize_from_params(params.get('initial_conditions', {}))
    
    # Print initial state
    print("\nInitial state:")
    print(f"  ne: min={state.electrons.n.x.array.min():.3e}, max={state.electrons.n.x.array.max():.3e} m⁻³")
    print(f"  Te: min={state.electrons.T.min():.3e}, max={state.electrons.T.max():.3e} eV")
    if 'Hi' in state.species:
        print(f"  nHi: min={state.species['Hi'].n.x.array.min():.3e}, max={state.species['Hi'].n.x.array.max():.3e} m⁻³")
    if 'HeII' in state.species:
        print(f"  nHeII: min={state.species['HeII'].n.x.array.min():.3e}, max={state.species['HeII'].n.x.array.max():.3e} m⁻³")
    
    # Create boundary handler
    bc_handler = BoundaryConditions(mesh, state.V)
    
    # Create transport
    transport = TransportManager(state.V)
    transport.initialize_from_params(
        params.get('transport', {}),
        [s.name for s in state.all_species]
    )
    
    # Create solver
    solver = BDF2Solver(state, transport, bc_handler, params.get('time_step', {}))
    
    # Set Dirichlet BC values for neutrals (fixed edge density)
    if 'nHeI_bc' in params:
        solver.set_dirichlet_value('HeI', params['nHeI_bc'])
        print(f"  HeI Dirichlet BC: {params['nHeI_bc']:.3e} m⁻³")
    if 'nH2_bc' in params:
        solver.set_dirichlet_value('H2', params['nH2_bc'])
        print(f"  H2 Dirichlet BC: {params['nH2_bc']:.3e} m⁻³")
    
    print(f"\nSolver initialized with dt={solver.dt:.3e} s")
    
    # Run a few debug steps
    print("\n" + "="*70)
    print("Running 20 debug steps...")
    print("="*70)
    
    for i in range(20):
        try:
            solver.step(debug=True)
        except Exception as e:
            print(f"\n!!! Exception at step {i}: {e}")
            import traceback
            traceback.print_exc()
            break
    
    print("\n" + "="*70)
    print("Final state:")
    print(f"  ne: min={state.electrons.n.x.array.min():.3e}, max={state.electrons.n.x.array.max():.3e} m⁻³")
    print(f"  Te: min={state.electrons.T.min():.3e}, max={state.electrons.T.max():.3e} eV")


if __name__ == "__main__":
    main()
