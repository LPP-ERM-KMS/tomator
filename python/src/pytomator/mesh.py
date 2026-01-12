"""
Mesh creation for 1D cylindrical plasma transport.

Supports loading radial positions from JSON file or creating uniform grids.
"""

import json
import numpy as np
from mpi4py import MPI
from dolfinx import mesh as dmesh
from dolfinx.mesh import Mesh


def create_mesh_from_json(filename: str, comm=MPI.COMM_WORLD) -> tuple[Mesh, np.ndarray]:
    """
    Create 1D interval mesh from JSON file containing radial positions.
    
    Expected JSON format:
    {
        "radial_positions": [r0, r1, r2, ..., rN]
    }
    
    Parameters
    ----------
    filename : str
        Path to JSON file with radial positions.
    comm : MPI communicator
        MPI communicator (default: COMM_WORLD).
        
    Returns
    -------
    mesh : dolfinx.mesh.Mesh
        1D interval mesh.
    radial_positions : np.ndarray
        Array of radial positions (for reference/output).
    """
    with open(filename, 'r') as f:
        data = json.load(f)
    
    radial_positions = np.array(data["radial_positions"], dtype=np.float64)
    
    # Validate monotonicity
    if not np.all(np.diff(radial_positions) > 0):
        raise ValueError("Radial positions must be strictly increasing.")
    
    r_min = radial_positions[0]
    r_max = radial_positions[-1]
    num_cells = len(radial_positions) - 1
    
    # Create interval mesh
    # Note: dolfinx creates uniform mesh, we'll need to move nodes for non-uniform
    domain = dmesh.create_interval(comm, num_cells, [r_min, r_max])
    
    # Move mesh nodes to match specified radial positions (for non-uniform grids)
    _adjust_mesh_coordinates(domain, radial_positions)
    
    return domain, radial_positions


def create_uniform_mesh(
    r_min: float, 
    r_max: float, 
    num_cells: int, 
    comm=MPI.COMM_WORLD
) -> tuple[Mesh, np.ndarray]:
    """
    Create uniform 1D interval mesh.
    
    Parameters
    ----------
    r_min : float
        Minimum radial position (HFS boundary).
    r_max : float
        Maximum radial position (LFS boundary).
    num_cells : int
        Number of mesh cells (num_points = num_cells + 1).
    comm : MPI communicator
        MPI communicator (default: COMM_WORLD).
        
    Returns
    -------
    mesh : dolfinx.mesh.Mesh
        1D interval mesh.
    radial_positions : np.ndarray
        Array of radial positions.
    """
    domain = dmesh.create_interval(comm, num_cells, [r_min, r_max])
    radial_positions = np.linspace(r_min, r_max, num_cells + 1)
    
    return domain, radial_positions


def create_mesh_from_geometry(
    R: float,
    a: float,
    lHFS: float,
    lLFS: float,
    num_cells: int,
    comm=MPI.COMM_WORLD
) -> tuple[Mesh, np.ndarray]:
    """
    Create mesh from tokamak geometry parameters.
    
    The grid spans from R - a to R + a (major radius ± minor radius).
    The limiter positions lHFS and lLFS define where the SOL begins
    and parallel losses occur (handled separately by limiters function).
    
    Parameters
    ----------
    R : float
        Major radius [m].
    a : float
        Minor radius [m].
    lHFS : float
        Absolute radial position of HFS limiter [m] (for parallel loss calc).
    lLFS : float
        Absolute radial position of LFS limiter [m] (for parallel loss calc).
    num_cells : int
        Number of mesh cells.
    comm : MPI communicator
        MPI communicator.
        
    Returns
    -------
    mesh : dolfinx.mesh.Mesh
        1D interval mesh.
    radial_positions : np.ndarray
        Array of radial positions.
    """
    r_min = R - a
    r_max = R + a
    
    return create_uniform_mesh(r_min, r_max, num_cells, comm)


def _adjust_mesh_coordinates(domain: Mesh, target_positions: np.ndarray) -> None:
    """
    Adjust mesh node coordinates to match target positions (for non-uniform grids).
    
    Parameters
    ----------
    domain : dolfinx.mesh.Mesh
        Mesh to modify in-place.
    target_positions : np.ndarray
        Target radial positions for each node.
    """
    # Get mesh geometry
    geometry = domain.geometry
    coords = geometry.x
    
    # Sort coordinates and map to target positions
    num_points = len(target_positions)
    
    if coords.shape[0] != num_points:
        raise ValueError(
            f"Mesh has {coords.shape[0]} points but {num_points} target positions provided."
        )
    
    # Get the sorting indices based on current x-coordinates
    sort_idx = np.argsort(coords[:, 0])
    
    # Assign target positions to sorted nodes
    for i, idx in enumerate(sort_idx):
        coords[idx, 0] = target_positions[i]


def get_boundary_facets(domain: Mesh) -> tuple[np.ndarray, np.ndarray]:
    """
    Get facet indices for HFS (left) and LFS (right) boundaries.
    
    Parameters
    ----------
    domain : dolfinx.mesh.Mesh
        1D interval mesh.
        
    Returns
    -------
    hfs_facets : np.ndarray
        Facet indices at HFS (minimum r) boundary.
    lfs_facets : np.ndarray
        Facet indices at LFS (maximum r) boundary.
    """
    # Get mesh coordinates
    coords = domain.geometry.x
    r_min = coords[:, 0].min()
    r_max = coords[:, 0].max()
    
    tol = 1e-10 * (r_max - r_min)
    
    # Locate boundary facets (vertices in 1D)
    def hfs_boundary(x):
        return np.isclose(x[0], r_min, atol=tol)
    
    def lfs_boundary(x):
        return np.isclose(x[0], r_max, atol=tol)
    
    fdim = domain.topology.dim - 1  # Facet dimension (0 for 1D)
    
    hfs_facets = dmesh.locate_entities_boundary(domain, fdim, hfs_boundary)
    lfs_facets = dmesh.locate_entities_boundary(domain, fdim, lfs_boundary)
    
    return hfs_facets, lfs_facets


def get_radial_coordinate(domain: Mesh):
    """
    Get the radial coordinate as a UFL expression.
    
    Parameters
    ----------
    domain : dolfinx.mesh.Mesh
        1D interval mesh.
        
    Returns
    -------
    r : ufl.SpatialCoordinate
        Radial coordinate for use in weak forms.
    """
    import ufl
    x = ufl.SpatialCoordinate(domain)
    return x[0]  # First (and only) coordinate in 1D
