"""
Boundary condition handling for plasma transport equations.

Supports:
- Robin (decay length) BCs: dn/dr = -n/lambda at boundaries
- Dirichlet BCs: fixed values at boundaries (for H2, HeI neutrals)
"""

from typing import Callable, Optional, Tuple
import numpy as np
import ufl
from dolfinx import fem, default_scalar_type
from dolfinx.mesh import Mesh, locate_entities_boundary, meshtags
from dolfinx.fem import Function, FunctionSpace, Constant, dirichletbc, locate_dofs_topological


class BoundaryConditions:
    """
    Container for boundary conditions on a 1D domain.
    
    Handles both HFS (high-field side, r_min) and LFS (low-field side, r_max)
    boundaries with either Robin or Dirichlet conditions.
    
    Attributes
    ----------
    mesh : dolfinx.mesh.Mesh
        The computational mesh.
    V : fem.FunctionSpace
        Function space.
    hfs_facets : np.ndarray
        Facet indices at HFS boundary.
    lfs_facets : np.ndarray
        Facet indices at LFS boundary.
    """
    
    def __init__(self, mesh: Mesh, V: FunctionSpace):
        """
        Initialize boundary condition handler.
        
        Parameters
        ----------
        mesh : dolfinx.mesh.Mesh
            1D interval mesh.
        V : fem.FunctionSpace
            Function space for the problem.
        """
        self.mesh = mesh
        self.V = V
        self.fdim = mesh.topology.dim - 1  # Facet dimension (0 for 1D)
        
        # Get boundary coordinates
        coords = mesh.geometry.x
        self.r_min = coords[:, 0].min()
        self.r_max = coords[:, 0].max()
        tol = 1e-10 * (self.r_max - self.r_min)
        
        # Locate boundary facets
        self.hfs_facets = locate_entities_boundary(
            mesh, self.fdim, 
            lambda x: np.isclose(x[0], self.r_min, atol=tol)
        )
        self.lfs_facets = locate_entities_boundary(
            mesh, self.fdim,
            lambda x: np.isclose(x[0], self.r_max, atol=tol)
        )
        
        # Create boundary measure for Robin BCs
        self._setup_boundary_measure()
    
    def _setup_boundary_measure(self) -> None:
        """Set up boundary facet tags and measure for weak form integration."""
        # Combine all boundary facets with markers
        # HFS = marker 1, LFS = marker 2
        all_facets = np.concatenate([self.hfs_facets, self.lfs_facets])
        markers = np.concatenate([
            np.full_like(self.hfs_facets, 1, dtype=np.int32),
            np.full_like(self.lfs_facets, 2, dtype=np.int32)
        ])
        
        # Sort by facet index
        sort_idx = np.argsort(all_facets)
        all_facets = all_facets[sort_idx]
        markers = markers[sort_idx]
        
        # Create meshtags
        self.mesh.topology.create_connectivity(self.fdim, self.mesh.topology.dim)
        self.facet_tags = meshtags(self.mesh, self.fdim, all_facets, markers)
        
        # Create boundary measure
        self.ds = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.facet_tags)
    
    def get_dirichlet_bc(
        self, 
        value: float, 
        boundary: str = "both"
    ) -> list:
        """
        Create Dirichlet boundary condition(s).
        
        Parameters
        ----------
        value : float
            Fixed value at boundary.
        boundary : str
            Which boundary: "hfs", "lfs", or "both".
            
        Returns
        -------
        bcs : list
            List of DirichletBC objects.
        """
        bcs = []
        
        if boundary in ("hfs", "both"):
            hfs_dofs = locate_dofs_topological(self.V, self.fdim, self.hfs_facets)
            bc_hfs = dirichletbc(
                default_scalar_type(value), 
                hfs_dofs, 
                self.V
            )
            bcs.append(bc_hfs)
        
        if boundary in ("lfs", "both"):
            lfs_dofs = locate_dofs_topological(self.V, self.fdim, self.lfs_facets)
            bc_lfs = dirichletbc(
                default_scalar_type(value),
                lfs_dofs,
                self.V
            )
            bcs.append(bc_lfs)
        
        return bcs
    
    def get_dirichlet_bc_function(
        self,
        value_func: Function,
        boundary: str = "both"
    ) -> list:
        """
        Create Dirichlet BC from a Function (for spatially varying values).
        
        Parameters
        ----------
        value_func : fem.Function
            Function containing boundary values.
        boundary : str
            Which boundary: "hfs", "lfs", or "both".
            
        Returns
        -------
        bcs : list
            List of DirichletBC objects.
        """
        bcs = []
        
        if boundary in ("hfs", "both"):
            hfs_dofs = locate_dofs_topological(self.V, self.fdim, self.hfs_facets)
            bc_hfs = dirichletbc(value_func, hfs_dofs)
            bcs.append(bc_hfs)
        
        if boundary in ("lfs", "both"):
            lfs_dofs = locate_dofs_topological(self.V, self.fdim, self.lfs_facets)
            bc_lfs = dirichletbc(value_func, lfs_dofs)
            bcs.append(bc_lfs)
        
        return bcs


class DecayLengthBC:
    """
    Robin boundary condition with decay length: dn/dr = -n/lambda.
    
    This is the standard BC for ions and most neutrals, where particles
    are lost to limiters with a characteristic decay length.
    
    The Robin BC is incorporated into the weak form as a boundary integral:
    
    For the equation: dn/dt = D * d²n/dr² + ...
    
    Weak form boundary term: + integral(D/lambda * n * v) ds
    
    This adds a loss term at the boundary proportional to the local density.
    
    Attributes
    ----------
    lambda_hfs : float
        Decay length at HFS boundary [m].
    lambda_lfs : float
        Decay length at LFS boundary [m].
    """
    
    def __init__(
        self,
        lambda_hfs: float,
        lambda_lfs: float,
        bc_handler: BoundaryConditions
    ):
        """
        Initialize decay length BC.
        
        Parameters
        ----------
        lambda_hfs : float
            Decay length at HFS (high-field side) boundary [m].
        lambda_lfs : float
            Decay length at LFS (low-field side) boundary [m].
        bc_handler : BoundaryConditions
            Boundary condition handler with mesh and measures.
        """
        self.lambda_hfs = lambda_hfs
        self.lambda_lfs = lambda_lfs
        self.bc_handler = bc_handler
        
        # Create constants for use in UFL forms
        self.inv_lambda_hfs = Constant(bc_handler.mesh, default_scalar_type(1.0 / lambda_hfs))
        self.inv_lambda_lfs = Constant(bc_handler.mesh, default_scalar_type(1.0 / lambda_lfs))
    
    def update_decay_lengths(self, lambda_hfs: float, lambda_lfs: float) -> None:
        """
        Update decay length values (e.g., temperature-dependent).
        
        Parameters
        ----------
        lambda_hfs : float
            New HFS decay length [m].
        lambda_lfs : float
            New LFS decay length [m].
        """
        self.lambda_hfs = lambda_hfs
        self.lambda_lfs = lambda_lfs
        self.inv_lambda_hfs.value = 1.0 / lambda_hfs
        self.inv_lambda_lfs.value = 1.0 / lambda_lfs
    
    def get_weak_form_terms(
        self,
        n: ufl.Argument,  # Trial function
        v: ufl.Argument,  # Test function
        D: ufl.Coefficient,  # Diffusion coefficient
        r: ufl.SpatialCoordinate,  # Radial coordinate for cylindrical geometry
        c1: ufl.Coefficient = None  # BDF2 coefficient for scaling spatial terms
    ) -> ufl.Form:
        """
        Get boundary integral terms for weak form.
        
        For cylindrical coordinates, the boundary terms include the r factor:
        
        a_boundary = c1 * integral(r * D/lambda_hfs * n * v) ds(1)
                   + c1 * integral(r * D/lambda_lfs * n * v) ds(2)
        
        Note: c1 coefficient is needed for proper BDF2 formulation where
        spatial terms (including Robin BC) are scaled by c1 while the
        mass term is not.
        
        Parameters
        ----------
        n : ufl.Argument
            Trial function (density).
        v : ufl.Argument
            Test function.
        D : ufl.Coefficient
            Diffusion coefficient (can be Function or Constant).
        r : ufl.SpatialCoordinate
            Radial coordinate.
        c1 : ufl.Coefficient, optional
            BDF2 coefficient for scaling. If None, no scaling is applied.
            
        Returns
        -------
        boundary_form : ufl.Form
            Boundary integral contribution to bilinear form.
        """
        ds = self.bc_handler.ds
        
        # HFS boundary (marker 1) - note: sign is + because Robin BC adds to LHS
        a_hfs = r * D * self.inv_lambda_hfs * n * v * ds(1)
        
        # LFS boundary (marker 2)
        a_lfs = r * D * self.inv_lambda_lfs * n * v * ds(2)
        
        boundary_form = a_hfs + a_lfs
        
        # Apply c1 scaling if provided (for BDF2 formulation)
        if c1 is not None:
            boundary_form = c1 * boundary_form
        
        return boundary_form
    
    def get_rhs_terms(
        self,
        v: ufl.Argument,
        D: ufl.Coefficient,
        r: ufl.SpatialCoordinate,
        n_bc: Optional[float] = None
    ) -> Optional[ufl.Form]:
        """
        Get RHS boundary terms if using inhomogeneous Robin BC.
        
        For homogeneous Robin BC (n -> 0 at boundary), this returns None.
        For inhomogeneous Robin BC (n -> n_bc at boundary), returns the
        corresponding RHS term.
        
        Parameters
        ----------
        v : ufl.Argument
            Test function.
        D : ufl.Coefficient
            Diffusion coefficient.
        r : ufl.SpatialCoordinate
            Radial coordinate.
        n_bc : float, optional
            Target boundary value (None for homogeneous BC).
            
        Returns
        -------
        rhs_form : ufl.Form or None
            RHS boundary contribution, or None for homogeneous BC.
        """
        if n_bc is None or n_bc == 0.0:
            return None
        
        ds = self.bc_handler.ds
        n_bc_const = Constant(self.bc_handler.mesh, default_scalar_type(n_bc))
        
        # RHS terms for inhomogeneous Robin BC
        L_hfs = r * D * self.inv_lambda_hfs * n_bc_const * v * ds(1)
        L_lfs = r * D * self.inv_lambda_lfs * n_bc_const * v * ds(2)
        
        return L_hfs + L_lfs


def compute_decay_length(
    T: float,
    m: float,
    connection_length: float = 10.0,
    cs_factor: float = 1.0
) -> float:
    """
    Compute decay length from temperature and connection length.
    
    lambda = D / (cs * something) ≈ rho_s or connection_length based
    
    For now, use a simplified model:
    lambda = connection_length * sqrt(m_i / m_e) * (T_i / T_e)^0.5
    
    Parameters
    ----------
    T : float
        Species temperature [eV].
    m : float
        Species mass [amu].
    connection_length : float
        Parallel connection length to limiter [m].
    cs_factor : float
        Sound speed correction factor.
        
    Returns
    -------
    decay_length : float
        Decay length [m].
    """
    # Sound speed: cs = sqrt(2 * T / m)
    EV_TO_J = 1.60218e-19
    AMU_TO_KG = 1.66054e-27
    
    T_J = T * EV_TO_J
    m_kg = m * AMU_TO_KG
    
    cs = np.sqrt(2 * T_J / m_kg) if T > 0 else 1.0
    
    # Simple model: lambda ~ D / cs, with D ~ rho_s * cs
    # For now, just use a fraction of connection length
    decay_length = connection_length * 0.01  # Placeholder - adjust based on physics
    
    # Ensure reasonable bounds
    decay_length = max(decay_length, 1e-4)  # Minimum 0.1 mm
    decay_length = min(decay_length, connection_length)
    
    return decay_length
