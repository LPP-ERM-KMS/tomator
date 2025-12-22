#!/usr/bin/env python3
"""
Example script demonstrating Tomator-dolfinx usage.

This example sets up a simple H + H+ + He plasma simulation
in a 1D cylindrical geometry, similar to the C++ Tomator1D code.
"""

import numpy as np
from mpi4py import MPI

# Import tomator_dolfinx modules
from tomator_dolfinx import (
    create_uniform_mesh,
    PlasmaState,
    BoundaryConditions,
    TransportManager,
    BDF2Solver,
    OutputManager,
)


def main():
    """Run a simple plasma transport simulation."""
    
    print("=" * 60)
    print("Tomator-dolfinx Example Simulation")
    print("=" * 60)
    
    # ===================================================================
    # 1. Define geometry and create mesh
    # ===================================================================
    
    # Tokamak-like geometry (in meters)
    R = 0.88          # Major radius [m]
    lHFS = 0.28       # Distance to HFS limiter [m]
    lLFS = 0.28       # Distance to LFS limiter [m]
    
    r_min = R - lHFS  # HFS boundary
    r_max = R + lLFS  # LFS boundary
    num_cells = 100   # Number of mesh cells
    
    print(f"\nGeometry:")
    print(f"  Major radius R = {R} m")
    print(f"  Radial domain: [{r_min:.3f}, {r_max:.3f}] m")
    print(f"  Mesh cells: {num_cells}")
    
    mesh, radial_positions = create_uniform_mesh(r_min, r_max, num_cells)
    
    # ===================================================================
    # 2. Create plasma state with species
    # ===================================================================
    
    state = PlasmaState(mesh, degree=1)
    
    # Add helium species (HeI with Dirichlet BC, HeII/HeIII with Robin)
    state.add_helium_species()
    
    print(f"\nSpecies: {[s.name for s in state.all_species]}")
    
    # ===================================================================
    # 3. Set initial conditions
    # ===================================================================
    
    initial_params = {
        # Electron temperature [eV]
        'Te0': 5.0,
        
        # Atomic hydrogen
        'nH0': 1e16,     # [m^-3]
        'TH0': 0.026,    # Room temperature [eV]
        
        # Protons (H+)
        'nHi0': 1e17,    # [m^-3]
        'THi0': 5.0,     # [eV]
        
        # Helium
        'nHeI0': 1e15,   # [m^-3]
        'THeI0': 0.026,  # [eV]
        'nHeII0': 1e14,  # [m^-3]
        'THeII0': 5.0,   # [eV]
        'nHeIII0': 1e12, # [m^-3]
        'THeIII0': 5.0,  # [eV]
    }
    
    state.initialize_from_params(initial_params)
    
    # Print initial state
    print(f"\nInitial conditions:")
    print(f"  Te = {initial_params['Te0']} eV")
    print(f"  n_e = {state.electrons.n.x.array.mean():.2e} m^-3 (from quasi-neutrality)")
    print(f"  n_H+ = {state.species['Hi'].n.x.array.mean():.2e} m^-3")
    print(f"  n_H = {state.species['H'].n.x.array.mean():.2e} m^-3")
    
    # ===================================================================
    # 4. Set up boundary conditions
    # ===================================================================
    
    bc_handler = BoundaryConditions(mesh, state.V)
    
    # ===================================================================
    # 5. Set up transport coefficients
    # ===================================================================
    
    transport = TransportManager(state.V)
    
    # Use fixed diffusion and convection (D and V from input)
    D_fix = 0.1   # [m²/s]
    V_fix = 0.01  # [m/s]
    
    transport_params = {
        'Dfix': D_fix,
        'Vfix': V_fix,
        'diffusion_type': 'fixed',
        'convection_type': 'fixed',
    }
    
    transport.initialize_from_params(
        transport_params,
        [s.name for s in state.all_species]
    )
    
    print(f"\nTransport:")
    print(f"  D = {D_fix} m²/s (fixed)")
    print(f"  V = {V_fix} m/s (fixed)")
    
    # ===================================================================
    # 6. Set up solver
    # ===================================================================
    
    solver_params = {
        'dtinit': 1e-8,         # Initial time step [s]
        'dtmin': 1e-10,         # Minimum time step [s]
        'dtmax': 1e-5,          # Maximum time step [s]
        'accur': 0.05,          # Target accuracy
        'decay_length_hfs': 0.02,  # HFS decay length [m]
        'decay_length_lfs': 0.02,  # LFS decay length [m]
        
        # Physics flags
        'bH': True,
        'bHe': True,
        'bcx': True,
        'bADAS': True,
    }
    
    solver = BDF2Solver(state, transport, bc_handler, solver_params)
    
    # Set Dirichlet BC value for HeI (neutral helium at edge)
    solver.set_dirichlet_value('HeI', initial_params['nHeI0'])
    
    print(f"\nSolver:")
    print(f"  Time step: [{solver_params['dtmin']:.0e}, {solver_params['dtmax']:.0e}] s")
    print(f"  Target accuracy: {solver_params['accur']}")
    
    # ===================================================================
    # 7. Set up output
    # ===================================================================
    
    output_dir = "results_example"
    output_manager = OutputManager(
        output_dir=output_dir,
        save_interval=1e-5,
        radial_positions=radial_positions
    )
    
    print(f"\nOutput directory: {output_dir}")
    
    # ===================================================================
    # 8. Run simulation
    # ===================================================================
    
    t_end = 1e-4  # End time [s]
    
    print(f"\n{'=' * 60}")
    print(f"Running simulation until t = {t_end:.0e} s")
    print(f"{'=' * 60}\n")
    
    step = 0
    print_interval = 100  # Print every N steps
    
    while solver.t < t_end:
        dt = solver.step()
        step += 1
        
        # Save output
        output_manager.save(state, solver.t)
        
        # Print progress
        if step % print_interval == 0:
            ne_center = state.electrons.n.x.array[len(radial_positions)//2]
            Te_center = state.electrons.T[len(radial_positions)//2]
            print(f"Step {step:5d}: t = {solver.t:.3e} s, dt = {dt:.3e} s, "
                  f"n_e = {ne_center:.2e} m^-3, T_e = {Te_center:.2f} eV")
    
    # Save final state
    output_manager.save_final(state, solver.t)
    
    # ===================================================================
    # 9. Print final results
    # ===================================================================
    
    print(f"\n{'=' * 60}")
    print("Simulation completed!")
    print(f"{'=' * 60}")
    print(f"\nFinal state at t = {solver.t:.3e} s:")
    print(f"  Total steps: {step}")
    print(f"  n_e (center) = {state.electrons.n.x.array[num_cells//2]:.2e} m^-3")
    print(f"  T_e (center) = {state.electrons.T[num_cells//2]:.2f} eV")
    print(f"  n_H+ (center) = {state.species['Hi'].n.x.array[num_cells//2]:.2e} m^-3")
    print(f"  n_HeII (center) = {state.species['HeII'].n.x.array[num_cells//2]:.2e} m^-3")
    
    print(f"\nOutput files written to: {output_dir}/")
    
    return state


if __name__ == "__main__":
    main()
