"""
Tomator-dolfinx: 1D plasma transport solver using FEniCS/dolfinx.

This package solves coupled reaction-diffusion-advection equations
for plasma species (ions, neutrals, electrons) in cylindrical geometry.

Example usage:
    from tomator_dolfinx import run_simulation
    
    # Run from input file
    state = run_simulation("input.json", output_dir="results/")
    
    # Or build manually:
    from tomator_dolfinx import (
        create_uniform_mesh, PlasmaState, BDF2Solver,
        BoundaryConditions, TransportManager
    )
    
    mesh, r = create_uniform_mesh(0.6, 1.0, 100)
    state = PlasmaState(mesh)
    state.add_helium_species()
    state.initialize_from_params({"Te0": 10.0, "nHi0": 1e18})
    # ... set up solver and run
"""

from .mesh import (
    create_mesh_from_json, 
    create_uniform_mesh,
    create_mesh_from_geometry,
    get_boundary_facets,
    get_radial_coordinate
)
from .species import Species, PlasmaState
from .boundary import BoundaryConditions, DecayLengthBC
from .transport import TransportCoefficients, TransportManager, DiffusionModel, ConvectionModel
from .solver import BDF2Solver, TransportEquation, run_simulation
from .io.json_input import load_input_file, create_default_params
from .io.output import write_csv_output, OutputManager
from .reactions import ReactionRates, compute_collision_sources

__version__ = "0.1.0"
__all__ = [
    # Mesh
    "create_mesh_from_json",
    "create_uniform_mesh",
    "create_mesh_from_geometry",
    "get_boundary_facets",
    "get_radial_coordinate",
    # Species
    "Species",
    "PlasmaState",
    # Boundary conditions
    "BoundaryConditions",
    "DecayLengthBC",
    # Transport
    "TransportCoefficients",
    "TransportManager",
    "DiffusionModel",
    "ConvectionModel",
    # Solver
    "BDF2Solver",
    "TransportEquation",
    "run_simulation",
    # I/O
    "load_input_file",
    "create_default_params",
    "write_csv_output",
    "OutputManager",
    # Reactions
    "ReactionRates",
    "compute_collision_sources",
]
