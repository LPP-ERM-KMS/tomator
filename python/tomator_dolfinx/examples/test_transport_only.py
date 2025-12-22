#!/usr/bin/env python3
"""
Test transport only (no reactions) to verify FEM solver stability.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import numpy as np
from tomator_dolfinx.io.json_input import load_input_file
from tomator_dolfinx.mesh import create_mesh_from_geometry
from tomator_dolfinx.species import PlasmaState
from tomator_dolfinx.boundary import BoundaryConditions
from tomator_dolfinx.transport import TransportManager
from tomator_dolfinx.solver import TransportEquation


def main():
    # Load TCV input file
    input_file = Path(__file__).parent / "TCV5151X_fixneDV_noH2.json"
    params = load_input_file(str(input_file))
    
    # Create mesh
    mesh, radial_positions = create_mesh_from_geometry(
        R=params.get('R', 1.0),
        a=params.get('a', 0.2),
        lHFS=params.get('lHFS', 0.2),
        lLFS=params.get('lLFS', 0.2),
        num_cells=params.get('nmeshp', 100) - 1
    )
    
    print(f"Mesh: r_min={radial_positions[0]:.4f} m, r_max={radial_positions[-1]:.4f} m")
    
    # Create state
    state = PlasmaState(mesh)
    if params.get('bHe', False):
        state.add_helium_species()
    state.initialize_from_params(params.get('initial_conditions', {}))
    
    # Create boundary handler
    bc_handler = BoundaryConditions(mesh, state.V)
    
    # Create transport equation solver
    transport_eq = TransportEquation(
        state.V, bc_handler, 
        decay_length_hfs=params.get('decay_length_hfs', 0.06),
        decay_length_lfs=params.get('decay_length_lfs', 0.11)
    )
    
    # Transport coefficients
    from dolfinx.fem import Function
    D = Function(state.V)
    D.x.array[:] = 0.333  # m²/s
    V = Function(state.V)
    V.x.array[:] = 3.33  # m/s
    
    # Zero source
    source = Function(state.V)
    source.x.array[:] = 0.0
    
    # Test species: HeII
    species = state.species['HeII']
    
    print(f"\nInitial HeII: min={species.n.x.array.min():.3e}, max={species.n.x.array.max():.3e}")
    
    # Time step
    dt = 1e-9
    transport_eq.update_bdf2_coefficients(dt, 0.0)
    
    # Run 10 steps with pure transport (no source)
    for i in range(10):
        # Solve (backward Euler first, then BDF2)
        if i == 0:
            transport_eq.update_bdf2_coefficients(dt, 0.0)
        else:
            transport_eq.update_bdf2_coefficients(dt, dt)
        
        transport_eq.solve(D, V, species.n, species.n_prev, species.n_prev2, source)
        
        # Store previous
        species.n_prev2.x.array[:] = species.n_prev.x.array
        species.n_prev.x.array[:] = species.n.x.array
        
        print(f"Step {i+1}: HeII min={species.n.x.array.min():.3e}, max={species.n.x.array.max():.3e}")
    
    print("\nTransport-only test complete.")


if __name__ == "__main__":
    main()
