"""
Main solver for 1D cylindrical plasma transport using dolfinx.

Implements:
- Weak form assembly for reaction-diffusion-advection equations
- BDF2 adaptive time stepping
- Species transport with Robin/Dirichlet boundary conditions
"""

from datetime import datetime
from typing import Optional, Dict, Tuple, Callable
import numpy as np
import ufl
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import fem, default_scalar_type
from dolfinx.fem import Function, FunctionSpace, Constant
from dolfinx.fem.petsc import assemble_matrix, assemble_vector, apply_lifting, set_bc

from .species import Species, PlasmaState
from .boundary import BoundaryConditions, DecayLengthBC, compute_physics_decay_length
from .transport import (TransportCoefficients, TransportManager, 
                        compute_neutral_diffusion_from_state, compute_neutral_diffusion)
from .reactions.collisions import compute_sources_from_state
from .parallel import compute_limiter_losses, compute_bpol_losses


class TransportEquation:
    """
    Weak form for 1D cylindrical transport equation.
    
    Solves:
    ∂n/∂t = (1/r) ∂/∂r(r D ∂n/∂r) - (1/r) ∂/∂r(r V n) + S
    
    In weak form with test function v and cylindrical r-weighting:
    ∫ r (n - n_prev)/dt v dr + ∫ r D ∇n·∇v dr + ∫ r V ∇n·v dr 
    + boundary terms = ∫ r S v dr
    """
    
    def __init__(
        self,
        V: FunctionSpace,
        bc_handler: BoundaryConditions,
        decay_length_hfs: float = 0.01,
        decay_length_lfs: float = 0.01
    ):
        """
        Initialize transport equation solver.
        
        Parameters
        ----------
        V : FunctionSpace
            Function space for the solution.
        bc_handler : BoundaryConditions
            Boundary condition handler.
        decay_length_hfs : float
            Decay length at HFS boundary [m].
        decay_length_lfs : float
            Decay length at LFS boundary [m].
        """
        self.V = V
        self.mesh = V.mesh
        self.bc_handler = bc_handler
        
        # Create Robin BC handler
        self.robin_bc = DecayLengthBC(decay_length_hfs, decay_length_lfs, bc_handler)
        
        # Get radial coordinate
        x = ufl.SpatialCoordinate(self.mesh)
        self.r = x[0]
        
        # Trial and test functions
        self.n_trial = ufl.TrialFunction(V)
        self.v_test = ufl.TestFunction(V)
        
        # Time step (will be updated)
        self.dt = Constant(self.mesh, default_scalar_type(1e-6))
        
        # BDF2 coefficients (will be updated)
        self.c1 = Constant(self.mesh, default_scalar_type(1.0))  # Implicit weight
        self.c2 = Constant(self.mesh, default_scalar_type(1.0))  # n^k weight  
        self.c3 = Constant(self.mesh, default_scalar_type(0.0))  # n^{k-1} weight
    
    def assemble_forms(
        self,
        D: Function,
        V_conv: Function,
        n_prev: Function,
        n_prev2: Function,
        source: Function,
        use_robin_bc: bool = True,
        destruction_rate: Function = None
    ) -> Tuple[ufl.Form, ufl.Form]:
        """
        Assemble bilinear and linear forms for the transport equation.
        
        BDF2 discretization:
        (c1 * n^{k+1} - c2 * n^k + c3 * n^{k-1}) / dt = RHS
        
        For stiff source terms, we split: S = S^+ - k_d * n
        where k_d is the destruction rate coefficient [1/s].
        S^+ (creation) stays explicit on RHS, k_d * n goes to LHS implicitly.
        
        Parameters
        ----------
        D : Function
            Diffusion coefficient [m²/s].
        V_conv : Function
            Convection velocity [m/s].
        n_prev : Function
            Solution at previous time step.
        n_prev2 : Function
            Solution at two time steps back.
        source : Function
            Source term [m^-3/s] (creation terms only if destruction_rate provided).
        use_robin_bc : bool
            Use Robin BC (True) or prepare for Dirichlet (False).
        destruction_rate : Function, optional
            Implicit destruction rate coefficient k_d [1/s].
            If provided, adds k_d * n term to LHS for implicit stability.
            
        Returns
        -------
        a : ufl.Form
            Bilinear form (LHS).
        L : ufl.Form
            Linear form (RHS).
        """
        r = self.r
        n = self.n_trial
        v = self.v_test
        dt = self.dt
        c1, c2, c3 = self.c1, self.c2, self.c3
        
        # Volume measure
        dx = ufl.dx
        
        # CORRECTED BDF2 weak form (matching C++ Tomator):
        # The C++ formulation is: (Ts - c1*Ds*tstep + c1*Vs*tstep) * n = Ts*c2*n_prev - Ts*c3*n_prev2 + ...
        # This means:
        #   - Mass term (Ts) is NOT scaled by c1
        #   - Spatial terms (Ds, Vs) ARE scaled by c1
        #   - RHS history terms have c2/dt and c3/dt coefficients
        #
        # Weak form: (1/dt)*n*v + c1*(diffusion + convection) = (c2/dt)*n_prev - (c3/dt)*n_prev2 + c1*source
        
        # Time derivative: r * (1/dt) * n * v (NOT scaled by c1!)
        a_time = r * n / dt * v * dx
        
        # Diffusion: c1 * r * D * ∇n · ∇v (integration by parts of ∇·(r D ∇n))
        # SCALED by c1 to match C++ Tomator
        a_diff = c1 * r * D * ufl.inner(ufl.grad(n), ufl.grad(v)) * dx
        
        # Convection: weak form of -∇·(V n) = -V·∇n - n·∇V
        # For constant V, this is -V·∇n
        # On LHS this should be NEGATIVE: -∫ r V ∂n/∂r v dr
        # SCALED by c1 to match C++ Tomator
        a_conv = -c1 * r * V_conv * ufl.grad(n)[0] * v * dx
        
        # Combine bilinear form
        a = a_time + a_diff + a_conv
        
        # Add implicit destruction term: k_d * n on LHS for stability
        # This treats -k_d * n implicitly instead of explicitly
        # SCALED by c1 to match spatial operator scaling
        if destruction_rate is not None:
            a_destruction = c1 * r * destruction_rate * n * v * dx
            a = a + a_destruction
        
        # Add Robin BC boundary terms if requested
        # SCALED by c1 to match spatial operator scaling
        if use_robin_bc:
            a_boundary = self.robin_bc.get_weak_form_terms(n, v, D, r, c1)
            a = a + a_boundary
        
        # Linear form (RHS): BDF2 time terms + source
        # BDF2: (c2 * n^k - c3 * n^{k-1}) / dt
        L_time = r * (c2 * n_prev - c3 * n_prev2) / dt * v * dx
        
        # Source term - SCALED by c1 to match C++ Tomator:
        # RHS has: coef1 * Ts * Snion * tstep, which in weak form is c1 * r * source * v
        L_source = c1 * r * source * v * dx
        
        # Combine linear form
        L = L_time + L_source
        
        return a, L
    
    def solve(
        self,
        D: Function,
        V_conv: Function,
        n: Function,
        n_prev: Function,
        n_prev2: Function,
        source: Function,
        bcs: list = None,
        use_robin_bc: bool = True,
        destruction_rate: Function = None
    ) -> None:
        """
        Solve the transport equation for one time step.
        
        Parameters
        ----------
        D : Function
            Diffusion coefficient.
        V_conv : Function
            Convection velocity.
        n : Function
            Solution function (will be updated in-place).
        n_prev : Function
            Solution at previous time step.
        n_prev2 : Function
            Solution at two time steps back.
        source : Function
            Source term (creation only if destruction_rate provided).
        bcs : list, optional
            List of Dirichlet BCs (for species that use Dirichlet instead of Robin).
        use_robin_bc : bool
            Whether to include Robin BC terms in weak form.
        destruction_rate : Function, optional
            Implicit destruction rate [1/s] to add to LHS for stability.
        """
        # Assemble forms
        a, L = self.assemble_forms(D, V_conv, n_prev, n_prev2, source, use_robin_bc,
                                   destruction_rate=destruction_rate)
        
        # Solve using direct PETSc assembly (compatible with dolfinx 0.10+)
        if bcs is None:
            bcs = []
        
        # Compile forms
        a_compiled = fem.form(a)
        L_compiled = fem.form(L)
        
        # Assemble matrix and vector
        A = assemble_matrix(a_compiled, bcs=bcs)
        A.assemble()
        
        b = assemble_vector(L_compiled)
        apply_lifting(b, [a_compiled], [bcs])
        b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
        set_bc(b, bcs)
        
        # Create solution vector
        n_new = fem.Function(self.V)
        
        # Solve with direct LU (fast for small 1D problems)
        ksp = PETSc.KSP().create(self.mesh.comm)
        ksp.setOperators(A)
        ksp.setType(PETSc.KSP.Type.PREONLY)
        ksp.getPC().setType(PETSc.PC.Type.LU)
        ksp.setFromOptions()
        ksp.solve(b, n_new.x.petsc_vec)
        n_new.x.scatter_forward()
        
        # Copy result to n
        n.x.array[:] = n_new.x.array
        
        # Clean up
        ksp.destroy()
        A.destroy()
        b.destroy()
    
    def update_dt(self, dt_new: float) -> None:
        """Update time step value."""
        self.dt.value = dt_new
    
    def update_bdf2_coefficients(self, dt_new: float, dt_old: float, step_count: int = 100) -> None:
        """
        Update BDF2 coefficients for variable time step.
        
        For BDF2 with variable dt:
        c1 = (1 + w) / (1 + 2w)
        c2 = (1 + w)² / (1 + 2w)  
        c3 = w² / (1 + 2w)
        
        where w = dt_new / dt_old
        
        Following C++ Tomator1D: use fixed BDF2 coefficients for first 10 steps
        (Nit <= 10) for startup stability.
        
        NOTE: With the CORRECTED BDF2 formulation (mass term NOT scaled by c1),
        the scheme is stable even from step 0. The fixed coefficient period is
        kept for consistency with C++ Tomator.
        
        Parameters
        ----------
        dt_new : float
            New time step.
        dt_old : float
            Previous time step.
        step_count : int
            Current step number (0-indexed). C++ uses Nit (1-indexed iteration counter).
        """
        # C++ uses Nit which starts at 1, so step_count < 10 corresponds to Nit <= 10
        if step_count < 10:
            # Fixed BDF2 coefficients for startup stability (matching C++ Nit <= 10)
            self.c1.value = 2.0 / 3.0
            self.c2.value = 4.0 / 3.0
            self.c3.value = 1.0 / 3.0
        elif dt_old > 0:
            # Variable BDF2 coefficients for adaptive time stepping
            w = dt_new / dt_old
            self.c1.value = (1 + w) / (1 + 2*w)
            self.c2.value = (1 + w)**2 / (1 + 2*w)
            self.c3.value = w**2 / (1 + 2*w)
        else:
            # Fallback to fixed BDF2 if dt_old is invalid
            self.c1.value = 2.0 / 3.0
            self.c2.value = 4.0 / 3.0
            self.c3.value = 1.0 / 3.0
        
        self.dt.value = dt_new


class BDF2Solver:
    """
    Adaptive BDF2 time integrator for multi-species plasma transport.
    
    Solves the coupled system of transport equations for all species
    using operator splitting (sequential solve for each species).
    """
    
    def __init__(
        self,
        state: PlasmaState,
        transport: TransportManager,
        bc_handler: BoundaryConditions,
        params: dict,
        radial_positions: np.ndarray = None
    ):
        """
        Initialize BDF2 solver.
        
        Parameters
        ----------
        state : PlasmaState
            Plasma state containing all species.
        transport : TransportManager
            Transport coefficient manager.
        bc_handler : BoundaryConditions
            Boundary condition handler.
        params : dict
            Simulation parameters including:
            - dtinit: Initial time step [s]
            - dtmin: Minimum time step [s]
            - dtmax: Maximum time step [s]
            - accur: Target accuracy (max relative change)
            - decay_length_hfs: HFS decay length [m]
            - decay_length_lfs: LFS decay length [m]
            - lHFS: HFS limiter position [m] (for parallel losses)
            - lLFS: LFS limiter position [m] (for parallel losses)
            - nlimiters: Number of poloidal limiters
        radial_positions : np.ndarray, optional
            Radial positions of mesh nodes [m]. If not provided, will
            be extracted from the mesh geometry.
        """
        self.state = state
        self.transport = transport
        self.bc_handler = bc_handler
        self.params = params
        
        # Store radial positions for limiter losses calculation
        if radial_positions is not None:
            self.radial_positions = radial_positions
        else:
            # Extract from mesh geometry
            self.radial_positions = state.mesh.geometry.x[:, 0].copy()
        
        # Time stepping parameters
        self.dt = params.get('dtinit', 1e-6)
        self.dt_prev = self.dt
        self.dt_min = params.get('dtmin', 1e-9)
        self.dt_max = params.get('dtmax', 1e-3)
        self.accuracy = params.get('accur', 0.05)
        
        # Decay lengths for Robin BC
        lambda_hfs = params.get('decay_length_hfs', 0.01)
        lambda_lfs = params.get('decay_length_lfs', 0.01)
        
        # Create transport equation solver
        self.transport_eq = TransportEquation(
            state.V, bc_handler, lambda_hfs, lambda_lfs
        )
        
        # Source term function
        self.source = Function(state.V, name="source")
        
        # Store Dirichlet BC values for species that need them
        self.dirichlet_values = {}
        
        # Current time
        self.t = 0.0
        
        # Step counter
        self.step_count = 0
        
        # Initialize coupled power if RF power is specified
        self.coupled_power = None
        self._pecabs = 0.0
        Prf = params.get('Prf', 0.0)
        if Prf > 0:
            from .coupledpower import CoupledPower
            self.coupled_power = CoupledPower(params)
    
    def set_dirichlet_value(self, species_name: str, value: float) -> None:
        """
        Set Dirichlet boundary value for a species.
        
        Parameters
        ----------
        species_name : str
            Name of species.
        value : float
            Boundary value.
        """
        self.dirichlet_values[species_name] = value

    def _update_robin_decay_lengths(
        self,
        species: Species,
        D: Function,
        m_amu: float = 1.0
    ) -> None:
        """
        Update Robin BC decay lengths based on local physics (C++ approach).
        
        Computes decay length from:
            λ = 2 * D / (vth * (1 - R))
        
        where:
            vth = sqrt(kB * T / m) is the thermal velocity
            R = reflection coefficient
        
        This is the physics-based approach used in C++ Tomator1D, where the
        decay length adapts to local temperature and diffusion coefficient
        rather than being a fixed geometric fraction.
        
        Parameters
        ----------
        species : Species
            Species to compute decay length for.
        D : Function
            Diffusion coefficient function.
        m_amu : float
            Particle mass in atomic mass units.
        """
        # Get reflection coefficient from params (default 0.5 for density)
        RH = self.params.get('RH', 0.5)
        
        coords = self.state.mesh.geometry.x[:, 0]
        D_arr = D.x.array
        
        # Get temperature at boundaries from energy/density
        E_arr = species.E.x.array
        n_arr = species.n.x.array
        
        sorted_idx = np.argsort(coords)
        
        # HFS boundary (smallest r)
        idx_hfs = sorted_idx[0]
        if n_arr[idx_hfs] > 1e-30:
            T_hfs = E_arr[idx_hfs] / (1.5 * n_arr[idx_hfs])
        else:
            # Use interior temperature as fallback
            idx_interior = sorted_idx[1]
            if n_arr[idx_interior] > 1e-30:
                T_hfs = E_arr[idx_interior] / (1.5 * n_arr[idx_interior])
            else:
                T_hfs = self.params.get('Ta0', 0.026)
        
        # LFS boundary (largest r)
        idx_lfs = sorted_idx[-1]
        if n_arr[idx_lfs] > 1e-30:
            T_lfs = E_arr[idx_lfs] / (1.5 * n_arr[idx_lfs])
        else:
            # Use interior temperature as fallback
            idx_interior = sorted_idx[-2]
            if n_arr[idx_interior] > 1e-30:
                T_lfs = E_arr[idx_interior] / (1.5 * n_arr[idx_interior])
            else:
                T_lfs = self.params.get('Ta0', 0.026)
        
        # Compute physics-based decay lengths
        lambda_hfs = compute_physics_decay_length(
            D=D_arr[idx_hfs],
            T=T_hfs,
            m_amu=m_amu,
            R=RH
        )
        
        lambda_lfs = compute_physics_decay_length(
            D=D_arr[idx_lfs],
            T=T_lfs,
            m_amu=m_amu,
            R=RH
        )
        
        # Update the Robin BC constants
        self.transport_eq.robin_bc.update_decay_lengths(lambda_hfs, lambda_lfs)

    def _update_robin_energy_decay_lengths(
        self,
        species: Species,
        D: Function,
        m_amu: float = 1.0
    ) -> None:
        """
        Update Robin BC decay lengths for energy equation (C++ approach).
        
        The C++ energy BC for H is:
            dE/dr = ±1/2 * gEdn * vth/D * n * 3/2 * T * (1 - REH)
        
        This uses REH (energy reflection) instead of RH, and includes gEdn factor.
        The energy decay length relates to density decay length by:
            λ_E = λ_n * (1 - RH) / (gEdn * (1 - REH))
        
        Parameters
        ----------
        species : Species
            Species to compute decay length for.
        D : Function
            Diffusion coefficient function.
        m_amu : float
            Particle mass in atomic mass units.
        """
        # Get parameters
        RH = self.params.get('RH', 0.5)
        REH = self.params.get('REH', 0.9)
        gEdn = self.params.get('gEdn', 5/3)
        
        coords = self.state.mesh.geometry.x[:, 0]
        D_arr = D.x.array
        
        # Get temperature at boundaries from energy/density
        E_arr = species.E.x.array
        n_arr = species.n.x.array
        
        sorted_idx = np.argsort(coords)
        
        # HFS boundary (smallest r)
        idx_hfs = sorted_idx[0]
        if n_arr[idx_hfs] > 1e-30:
            T_hfs = E_arr[idx_hfs] / (1.5 * n_arr[idx_hfs])
        else:
            idx_interior = sorted_idx[1]
            if n_arr[idx_interior] > 1e-30:
                T_hfs = E_arr[idx_interior] / (1.5 * n_arr[idx_interior])
            else:
                T_hfs = self.params.get('Ta0', 0.026)
        
        # LFS boundary (largest r)
        idx_lfs = sorted_idx[-1]
        if n_arr[idx_lfs] > 1e-30:
            T_lfs = E_arr[idx_lfs] / (1.5 * n_arr[idx_lfs])
        else:
            idx_interior = sorted_idx[-2]
            if n_arr[idx_interior] > 1e-30:
                T_lfs = E_arr[idx_interior] / (1.5 * n_arr[idx_interior])
            else:
                T_lfs = self.params.get('Ta0', 0.026)
        
        # Compute density-based decay lengths first
        lambda_n_hfs = compute_physics_decay_length(
            D=D_arr[idx_hfs],
            T=T_hfs,
            m_amu=m_amu,
            R=RH
        )
        
        lambda_n_lfs = compute_physics_decay_length(
            D=D_arr[idx_lfs],
            T=T_lfs,
            m_amu=m_amu,
            R=RH
        )
        
        # Update Robin BC with energy-specific decay lengths
        self.transport_eq.robin_bc.update_energy_decay_lengths(
            lambda_n_hfs, lambda_n_lfs, gEdn, REH, RH
        )

    def _compute_neutral_diffusion(
        self,
        species_name: str,
        D: Function
    ) -> None:
        """
        Compute physics-based diffusion coefficient for neutral species.
        
        Neutrals (H, H2, HeI) always use physics-based diffusion:
            D = (1/3) * vth / (1/mfp + 1/(a/2))
        
        This is the C++ Tomator1D approach where neutral diffusion is
        determined by collision physics.
        
        Parameters
        ----------
        species_name : str
            Neutral species name ('H', 'H2', or 'HeI').
        D : Function
            Diffusion coefficient function to update in-place.
        """
        # Only compute for neutrals
        if species_name not in ('H', 'H2', 'HeI'):
            return
        
        # Get geometry
        a_minor = self.params.get('a', 0.1)  # Minor radius [m]
        
        # Get nu_collision if available (computed in step())
        nu_collision = getattr(self, '_nu_collision', None)
        
        # Compute physics-based diffusion with self-collisions for neutrals
        D_physics = compute_neutral_diffusion_from_state(
            species_name=species_name,
            state=self.state,
            a_minor=a_minor,
            include_self_collision=True,
            nu_collision=nu_collision
        )
        
        # Update diffusion function
        D.x.array[:] = D_physics

    def _apply_energy_bc_with_flux(
        self, 
        E: Function, 
        n: Function,
        D: Function, 
        n_bc: float, 
        Ta0: float,
        REH: float = 0.9
    ) -> None:
        """
        Apply energy BC based on flux direction at boundaries.
        
        The flux direction determines the energy boundary condition:
        - Inward flux: incoming particles at Ta0 + reflected particles at REH * T_interior
        - Outward flux: particles leave with their local/interior temperature
        
        Flux: Γ = -D * dn/dr
        At HFS (r_min): Γ > 0 means flux toward +r = INTO domain
        At LFS (r_max): Γ < 0 means flux toward -r = INTO domain
        
        Parameters
        ----------
        E : Function
            Energy density function to modify at boundaries.
        n : Function
            Density function (after transport solve).
        D : Function
            Diffusion coefficient function.
        n_bc : float
            Dirichlet boundary value for density.
        Ta0 : float
            Ambient/wall temperature [eV].
        REH : float
            Energy reflection coefficient (0 to 1). Default 0.9.
        """
        coords = self.state.mesh.geometry.x[:, 0]  # Radial coordinates
        E_arr = E.x.array
        n_arr = n.x.array
        D_arr = D.x.array
        
        # Sort by radial position to find neighbors
        sorted_idx = np.argsort(coords)
        r_sorted = coords[sorted_idx]
        
        # Energy at ambient temperature for one particle: E = 3/2 * T
        E_ambient = 1.5 * Ta0
        
        # === HFS boundary (smallest r) ===
        idx_hfs = sorted_idx[0]           # Boundary DOF index
        idx_interior_hfs = sorted_idx[1]  # First interior point
        
        dr_hfs = r_sorted[1] - r_sorted[0]
        dn_dr_hfs = (n_arr[idx_interior_hfs] - n_arr[idx_hfs]) / dr_hfs
        Gamma_hfs = -D_arr[idx_hfs] * dn_dr_hfs  # Positive = toward +r = INTO domain
        
        n_interior_hfs = n_arr[idx_interior_hfs]
        E_interior_hfs = E_arr[idx_interior_hfs]
        
        if Gamma_hfs > 0:  # Inward flux at HFS
            # Incoming particles at Ta0, reflected particles at T_interior
            if n_interior_hfs > 1e-30:
                T_interior = E_interior_hfs / (1.5 * n_interior_hfs)
                T_effective = (1.0 - REH) * Ta0 + REH * T_interior
                E_arr[idx_hfs] = n_bc * 1.5 * T_effective
            else:
                E_arr[idx_hfs] = n_bc * E_ambient
        else:  # Outward flux at HFS
            # Particles leave with their interior temperature
            if n_interior_hfs > 1e-30:
                E_per_particle = E_interior_hfs / n_interior_hfs
                E_arr[idx_hfs] = n_bc * E_per_particle
            else:
                E_arr[idx_hfs] = n_bc * E_ambient
        
        # === LFS boundary (largest r) ===
        idx_lfs = sorted_idx[-1]           # Boundary DOF index
        idx_interior_lfs = sorted_idx[-2]  # First interior point
        
        dr_lfs = r_sorted[-1] - r_sorted[-2]
        dn_dr_lfs = (n_arr[idx_lfs] - n_arr[idx_interior_lfs]) / dr_lfs
        Gamma_lfs = -D_arr[idx_lfs] * dn_dr_lfs  # Positive = toward +r = OUT of domain
        
        n_interior_lfs = n_arr[idx_interior_lfs]
        E_interior_lfs = E_arr[idx_interior_lfs]
        
        if Gamma_lfs < 0:  # Inward flux at LFS (negative = toward -r = INTO domain)
            # Incoming particles at Ta0, reflected particles at T_interior
            if n_interior_lfs > 1e-30:
                T_interior = E_interior_lfs / (1.5 * n_interior_lfs)
                T_effective = (1.0 - REH) * Ta0 + REH * T_interior
                E_arr[idx_lfs] = n_bc * 1.5 * T_effective
            else:
                E_arr[idx_lfs] = n_bc * E_ambient
        else:  # Outward flux at LFS
            # Particles leave with their interior temperature
            if n_interior_lfs > 1e-30:
                E_per_particle = E_interior_lfs / n_interior_lfs
                E_arr[idx_lfs] = n_bc * E_per_particle
            else:
                E_arr[idx_lfs] = n_bc * E_ambient
    
    def step(self, debug: bool = False) -> float:
        """
        Perform one time step.
        
        Returns
        -------
        dt : float
            The time step that was taken.
        """
        # Debug helper
        def _dbg(msg, arrays=None):
            if not debug:
                return
            print(f"  [DBG] {msg}")
            if arrays:
                for name, arr in arrays.items():
                    arr_np = arr if isinstance(arr, np.ndarray) else arr.x.array
                    print(f"        {name}: min={arr_np.min():.3e}, max={arr_np.max():.3e}")
        
        if debug:
            print(f"\n=== STEP {self.step_count}, t={self.t:.3e}, dt={self.dt:.3e} ===")
            _dbg("Initial state:", {
                'ne': self.state.electrons.n,
                'nHi': self.state.species['Hi'].n if 'Hi' in self.state.species else None,
                'nHeI': self.state.species['HeI'].n if 'HeI' in self.state.species else None,
                'nHeII': self.state.species['HeII'].n if 'HeII' in self.state.species else None,
            })
        
        # 1. Update BDF2 coefficients (pass step_count for startup stability)
        self.transport_eq.update_bdf2_coefficients(self.dt, self.dt_prev, self.step_count)
        
        # 2. Compute collision source terms (dn, dE, nu)
        dn_sources, dE_sources, nu_collision = compute_sources_from_state(self.state, self.params)
        
        # 2b. Compute parallel transport losses to limiters (SOL region)
        lHFS = self.params.get('lHFS', None)
        lLFS = self.params.get('lLFS', None)
        nlimiters = self.params.get('nlimiters', 6)
        
        if lHFS is not None and lLFS is not None:
            # Extract densities and temperatures for limiter loss calculation
            ne = self.state.electrons.n.x.array.copy()
            Te = self.state.electrons.T.copy()
            nHi = self.state.species['Hi'].n.x.array.copy() if 'Hi' in self.state.species else np.zeros_like(ne)
            THi = self.state.species['Hi'].T.copy() if 'Hi' in self.state.species else Te.copy()
            
            # Optional molecular species
            nH2i = self.state.species['H2i'].n.x.array.copy() if 'H2i' in self.state.species else None
            TH2i = self.state.species['H2i'].T.copy() if 'H2i' in self.state.species else None
            nH3i = self.state.species['H3i'].n.x.array.copy() if 'H3i' in self.state.species else None
            TH3i = self.state.species['H3i'].T.copy() if 'H3i' in self.state.species else None
            
            # Optional helium species
            nHeII = self.state.species['HeII'].n.x.array.copy() if 'HeII' in self.state.species else None
            THeII = self.state.species['HeII'].T.copy() if 'HeII' in self.state.species else None
            nHeIII = self.state.species['HeIII'].n.x.array.copy() if 'HeIII' in self.state.species else None
            THeIII = self.state.species['HeIII'].T.copy() if 'HeIII' in self.state.species else None
            
            dn_limiter, dE_limiter = compute_limiter_losses(
                self.radial_positions,
                lHFS, lLFS, nlimiters,
                ne, Te, nHi, THi,
                nH2i, TH2i, nH3i, TH3i,
                nHeII, THeII, nHeIII, THeIII,
                Ta0=self.params.get('Ta0', 0.026)
            )
            
            # Add limiter losses to collision sources
            for species_name, dn_lim in dn_limiter.items():
                if species_name in dn_sources:
                    dn_sources[species_name] += dn_lim
                else:
                    dn_sources[species_name] = dn_lim
                    
            for species_name, dE_lim in dE_limiter.items():
                if species_name in dE_sources:
                    dE_sources[species_name] += dE_lim
                else:
                    dE_sources[species_name] = dE_lim
        
        # 2c. Compute vertical diffusion losses (bpol_function)
        Bt = self.params.get('Bt', None)
        Bv = self.params.get('Bv', None)
        b_vert = self.params.get('b', None)  # Vertical extent [cm]
        R0 = self.params.get('R', None)  # Major radius [cm]
        
        if Bt is not None and Bv is not None and b_vert is not None and R0 is not None:
            # Compute local toroidal field: Br = Bt * R0 / R [T]
            # R0 and radial_positions are in different units - need to be consistent
            # radial_positions is in [m], R0 is in [cm]
            R_cm = self.radial_positions * 100.0  # Convert m to cm
            Br = Bt * R0 / R_cm  # [T], 1/R scaling
            
            # Extract densities, temperatures, energies, and collision frequencies
            ne = self.state.electrons.n.x.array.copy()
            Te = self.state.electrons.T.copy()
            nHi = self.state.species['Hi'].n.x.array.copy() if 'Hi' in self.state.species else np.zeros_like(ne)
            THi = self.state.species['Hi'].T.copy() if 'Hi' in self.state.species else Te.copy()
            EHi = 1.5 * THi * nHi  # Energy density [eV·m^-3]
            nuHi = nu_collision.get('Hi', np.ones_like(ne) * 5e3)  # Use computed nu or minimum
            
            # Optional molecular species
            nH2i = self.state.species['H2i'].n.x.array.copy() if 'H2i' in self.state.species else None
            TH2i = self.state.species['H2i'].T.copy() if 'H2i' in self.state.species else None
            EH2i = 1.5 * TH2i * nH2i if nH2i is not None else None
            nuH2i = nu_collision.get('H2i', None)
            
            nH3i = self.state.species['H3i'].n.x.array.copy() if 'H3i' in self.state.species else None
            TH3i = self.state.species['H3i'].T.copy() if 'H3i' in self.state.species else None
            EH3i = 1.5 * TH3i * nH3i if nH3i is not None else None
            nuH3i = nu_collision.get('H3i', None)
            
            # Optional helium species
            nHeII = self.state.species['HeII'].n.x.array.copy() if 'HeII' in self.state.species else None
            THeII = self.state.species['HeII'].T.copy() if 'HeII' in self.state.species else None
            EHeII = 1.5 * THeII * nHeII if nHeII is not None else None
            nuHeII = nu_collision.get('HeII', None)
            
            nHeIII = self.state.species['HeIII'].n.x.array.copy() if 'HeIII' in self.state.species else None
            THeIII = self.state.species['HeIII'].T.copy() if 'HeIII' in self.state.species else None
            EHeIII = 1.5 * THeIII * nHeIII if nHeIII is not None else None
            nuHeIII = nu_collision.get('HeIII', None)
            
            # Get diffusion parameters
            Dfsave = self.params.get('Dfsave', 1.0)
            Dfix = self.params.get('Dfix', None) if self.params.get('bDfix', False) else None
            gEd = self.params.get('gEd', 1.0)
            
            dn_bpol, dE_bpol = compute_bpol_losses(
                Br=Br,
                Bv=Bv,
                b=b_vert,
                ne=ne, Te=Te,
                nHi=nHi, THi=THi, EHi=EHi, nuHi=nuHi,
                nH2i=nH2i, TH2i=TH2i, EH2i=EH2i, nuH2i=nuH2i,
                nH3i=nH3i, TH3i=TH3i, EH3i=EH3i, nuH3i=nuH3i,
                nHeII=nHeII, THeII=THeII, EHeII=EHeII, nuHeII=nuHeII,
                nHeIII=nHeIII, THeIII=THeIII, EHeIII=EHeIII, nuHeIII=nuHeIII,
                Dfsave=Dfsave,
                Dfix=Dfix,
                gEd=gEd,
                Ta0=self.params.get('Ta0', 0.026)
            )
            
            # Add bpol losses to collision sources
            for species_name, dn_bp in dn_bpol.items():
                if species_name in dn_sources:
                    dn_sources[species_name] += dn_bp
                else:
                    dn_sources[species_name] = dn_bp
                    
            for species_name, dE_bp in dE_bpol.items():
                if species_name in dE_sources:
                    dE_sources[species_name] += dE_bp
                else:
                    dE_sources[species_name] = dE_bp
        
        # 2d. Compute RF power deposition (coupled power)
        if hasattr(self, 'coupled_power') and self.coupled_power is not None:
            power_result = self.coupled_power.compute_power(
                R=self.radial_positions,
                ne=self.state.electrons.n.x.array.copy(),
                Te=self.state.electrons.T.copy(),
                dt=self.dt,
                t=self.t,
                nue=nu_collision.get('e', None)
            )
            
            # Add power to electron energy source
            # PRFe is in [eV/m³/s] (already SI since R and b are in meters)
            if 'e' in dE_sources:
                dE_sources['e'] += power_result.PRFe
            else:
                dE_sources['e'] = power_result.PRFe.copy()
            
            # Store absorbed power fraction for diagnostics
            self._pecabs = power_result.pecabs
        
        if debug:
            _dbg("Source terms (dn/dt):", {
                k: v for k, v in dn_sources.items() if v is not None
            })
        
        # Store nu_collision for use in transport calculations if needed
        self._nu_collision = nu_collision
        
        # 2e. Update transport coefficients (Bohm, classical, etc.)
        # For Bohm diffusion: D = Dfsave * T_eff / B
        # Compute radial magnetic field: Br = Bt * R0 / R
        Bt = self.params.get('Bt', None)
        R0 = self.params.get('R', None)  # Major radius [cm]
        
        if Bt is not None and R0 is not None:
            R_cm = self.radial_positions * 100.0  # Convert m to cm
            B_radial = Bt * R0 / R_cm  # [T], 1/R scaling
            self.transport.update_all(self.state, B_field=B_radial)
        else:
            # No magnetic field info - only FIXED diffusion will work
            self.transport.update_all(self.state, B_field=None)
        
        # 3. Solve transport for each species
        self._solve_ions(dn_sources, dE_sources)
        
        if debug:
            _dbg("After ion solve:", {
                'nHi': self.state.species['Hi'].n if 'Hi' in self.state.species else None,
                'nHeII': self.state.species['HeII'].n if 'HeII' in self.state.species else None,
            })
        
        self._solve_neutrals(dn_sources, dE_sources)
        
        if debug:
            _dbg("After neutral solve:", {
                'nH': self.state.species['H'].n if 'H' in self.state.species else None,
            })
        
        # 4. Apply density floors (CRITICAL for stability, matching C++)
        self._apply_density_floors()
        
        # 4b. Apply temperature clamps (CRITICAL for stability, matching C++)
        self._apply_temperature_clamps()
        
        # 5. Update electron density from quasi-neutrality
        self.state.compute_electron_density()
        
        if debug:
            _dbg("After QN ne update:", {
                'ne': self.state.electrons.n,
            })
        
        # 6. Solve electron energy equation
        self._solve_electron_energy(dE_sources)
        
        if debug:
            Te = self.state.electrons.T
            _dbg("After energy solve:", {
                'Te': Te,
            })
        
        # 7. Adapt time step
        dt_taken = self.dt
        self._adapt_timestep()
        
        # 8. Store previous solutions for BDF2
        self.state.store_all_previous()
        
        # 9. Update time and counter
        self.t += dt_taken
        self.step_count += 1
        self.dt_prev = dt_taken
        
        return dt_taken
    
    def _apply_density_floors(self) -> None:
        """
        Apply density floor clamping to all ion species.
        
        This matches C++ Tomator1D timeStep.cpp:
        - nevac = 1.0 cm⁻³ = 1e6 m⁻³ is the floor for all ion densities
        - When density is clamped to floor, energy is rescaled proportionally
        
        This is CRITICAL for numerical stability. Without it, densities can
        go negative or become extremely small, causing rate coefficient
        calculations to blow up.
        """
        # nevac from C++: 1.0 cm⁻³ = 1e6 m⁻³ in SI units
        nevac_si = self.params.get('nevac', 1.0) * 1e6  # Convert cm⁻³ to m⁻³
        
        # Apply floors to all ion species
        ion_names = ['Hi', 'H2i', 'H3i', 'HeII', 'HeIII']
        
        for name in ion_names:
            if name not in self.state.species:
                continue
                
            species = self.state.species[name]
            n_arr = species.n.x.array
            E_arr = species.E.x.array
            
            # Find where density is below floor
            below_floor = n_arr < nevac_si
            
            if np.any(below_floor):
                # Rescale energy proportionally when clamping density
                # E_new = E_old * nevac / n_old (to preserve temperature)
                scale_factor = np.where(below_floor & (n_arr > 1e-30),
                                       nevac_si / n_arr, 1.0)
                E_arr[:] = E_arr * scale_factor
                
                # Clamp density to floor
                n_arr[:] = np.maximum(n_arr, nevac_si)
        
        # Also apply floor to neutral species (optional, for stability)
        neutral_names = ['H', 'H2', 'HeI']
        for name in neutral_names:
            if name not in self.state.species:
                continue
            species = self.state.species[name]
            n_arr = species.n.x.array
            # Use smaller floor for neutrals since they can be very dilute
            n_arr[:] = np.maximum(n_arr, 1e-10)
    
    def _apply_temperature_clamps(self) -> None:
        """
        Clamp energy density to enforce temperature limits.
        
        This is CRITICAL for numerical stability when densities are low.
        The temperature T = E / (1.5 * n) can blow up if n is small.
        We enforce E <= 1.5 * n * T_max for all species.
        
        Temperature limits (matching collisions.py):
        - Heavy particles (ions, neutrals): T_max = 1000 eV
        - Electrons: T_max = 20000 eV
        """
        T_MIN = 0.026  # Room temperature
        T_MAX_HEAVY = 1000.0  # Max for ions/neutrals
        T_MAX_ELECTRON = 2e4  # Max for electrons
        
        # Apply to all heavy species
        all_names = ['Hi', 'H2i', 'H3i', 'HeII', 'HeIII', 'H', 'H2', 'HeI']
        for name in all_names:
            if name not in self.state.species:
                continue
            species = self.state.species[name]
            n_arr = species.n.x.array
            E_arr = species.E.x.array
            
            # E = 1.5 * n * T, so E_max = 1.5 * n * T_max
            E_min = 1.5 * n_arr * T_MIN
            E_max = 1.5 * n_arr * T_MAX_HEAVY
            E_arr[:] = np.clip(E_arr, E_min, E_max)
        
        # Apply to electrons
        if self.state.electrons is not None:
            n_arr = self.state.electrons.n.x.array
            E_arr = self.state.electrons.E.x.array
            E_min = 1.5 * n_arr * T_MIN
            E_max = 1.5 * n_arr * T_MAX_ELECTRON
            E_arr[:] = np.clip(E_arr, E_min, E_max)
    
    def _enforce_outflow_boundary(self, n: Function) -> None:
        """
        Enforce outflow-only boundary condition for ions.
        
        This matches C++ Tomator1D transport.cpp:
        - At HFS boundary: gradient must be >= 0 (density increasing outward from wall)
        - At LFS boundary: gradient must be <= 0 (density decreasing toward wall)
        
        If gradient has wrong sign (would cause inward flux), we clamp the
        boundary value to prevent artificial ion influx from outside domain.
        
        Parameters
        ----------
        n : Function
            Density function to enforce boundary constraints on.
        """
        coords = self.state.mesh.geometry.x[:, 0]  # Radial coordinates
        n_arr = n.x.array
        
        # Find boundary indices (assume sorted by radius)
        idx_hfs = np.argmin(coords)  # HFS = smallest r
        idx_lfs = np.argmax(coords)  # LFS = largest r
        
        # Get neighboring indices for gradient calculation
        # For HFS: check if density is decreasing into domain (bad)
        # For LFS: check if density is increasing into domain (bad)
        
        # Find second point from HFS (next larger r)
        if idx_hfs < len(coords) - 1:
            # Sort to find ordering
            sorted_idx = np.argsort(coords)
            idx_hfs_sorted = np.where(sorted_idx == idx_hfs)[0][0]
            if idx_hfs_sorted < len(sorted_idx) - 1:
                idx_next_hfs = sorted_idx[idx_hfs_sorted + 1]
                # Gradient at HFS: (n_next - n_hfs) / dr
                dr_hfs = coords[idx_next_hfs] - coords[idx_hfs]
                if dr_hfs > 0:
                    grad_hfs = (n_arr[idx_next_hfs] - n_arr[idx_hfs]) / dr_hfs
                    # If gradient < 0, density is higher at boundary than inside
                    # This would cause inward flux - clamp boundary to interior value
                    if grad_hfs < 0:
                        n_arr[idx_hfs] = n_arr[idx_next_hfs]
        
        # Find second point from LFS (next smaller r)
        if idx_lfs > 0:
            idx_lfs_sorted = np.where(sorted_idx == idx_lfs)[0][0]
            if idx_lfs_sorted > 0:
                idx_next_lfs = sorted_idx[idx_lfs_sorted - 1]
                # Gradient at LFS: (n_lfs - n_prev) / dr
                dr_lfs = coords[idx_lfs] - coords[idx_next_lfs]
                if dr_lfs > 0:
                    grad_lfs = (n_arr[idx_lfs] - n_arr[idx_next_lfs]) / dr_lfs
                    # If gradient > 0, density is higher at boundary than inside
                    # This would cause inward flux - clamp boundary to interior value
                    if grad_lfs > 0:
                        n_arr[idx_lfs] = n_arr[idx_next_lfs]
    
    def _solve_ions(self, dn_sources: dict, dE_sources: dict) -> None:
        """Solve transport equations for ion species."""
        # Map ions to their atomic masses (for decay length calculation)
        ion_mass = {'Hi': 1.0, 'H2i': 2.0, 'H3i': 3.0, 'HeII': 4.0, 'HeIII': 4.0}
        
        for species in self.state.ions:
            if not species.solve_density:
                continue
            
            name = species.name
            
            # Get transport coefficients (updated by transport.update_all in step())
            # D model depends on params: bDfix -> fixed, bDbohm -> Bohm, bDscaling -> classical
            coeff = self.transport.get(name)
            D = coeff.D
            V = coeff.V
            
            # Set source term
            if name in dn_sources:
                self.source.x.array[:] = dn_sources[name]
            else:
                self.source.x.array[:] = 0.0
            
            # Determine BC type
            use_robin = (species.bc_type == "robin")
            bcs = []
            
            if not use_robin and name in self.dirichlet_values:
                # Create Dirichlet BCs
                bcs = self.bc_handler.get_dirichlet_bc(
                    self.dirichlet_values[name], "both"
                )
            
            # For Robin BC species: update decay lengths from physics
            # λ = 2 * D / (vth * (1 - R))
            if use_robin:
                m_amu = ion_mass.get(name, 1.0)
                self._update_robin_decay_lengths(species, D, m_amu)
            
            # Solve density equation
            self.transport_eq.solve(
                D, V, species.n, species.n_prev, species.n_prev2,
                self.source, bcs=bcs, use_robin_bc=use_robin
            )
            
            # CRITICAL: Enforce outflow-only at boundaries for ions
            # This prevents artificial ion influx that causes numerical blowup
            self._enforce_outflow_boundary(species.n)
            
            # Solve energy equation (similar structure)
            # C++ uses: B = Ts - coef1 * gEd * Ds * tstep + coef1 * gEv * Vs * tstep
            # So for energy, scale D by gEd and V by gEv
            if species.solve_energy and name in dE_sources:
                self.source.x.array[:] = dE_sources[name]
                
                # Get energy transport scaling factors
                gEd = self.params.get('gEd', 5/3)
                gEv = self.params.get('gEv', 5/3)
                
                # Create scaled transport coefficients for energy
                D_energy = Function(self.state.V)
                V_energy = Function(self.state.V)
                D_energy.x.array[:] = gEd * D.x.array[:]
                V_energy.x.array[:] = gEv * V.x.array[:]
                
                # Use energy-specific decay lengths for Robin BC species
                if use_robin:
                    m_amu = ion_mass.get(name, 1.0)
                    self._update_robin_energy_decay_lengths(species, D, m_amu)
                
                self.transport_eq.solve(
                    D_energy, V_energy, species.E, species.E_prev, species.E_prev2,
                    self.source, bcs=[], use_robin_bc=use_robin
                )
                # Also enforce outflow for energy
                self._enforce_outflow_boundary(species.E)
    
    def _solve_neutrals(self, dn_sources: dict, dE_sources: dict) -> None:
        """Solve transport equations for neutral species."""
        # Get parameters for flux-based energy BC
        Ta0 = self.params.get('Ta0', 0.026)  # Ambient temperature [eV]
        REH = self.params.get('REH', 0.9)    # Energy reflection coefficient
        
        # Map species to their atomic masses (for decay length calculation)
        species_mass = {'H': 1.0, 'H2': 2.0, 'HeI': 4.0}
        
        for species in self.state.neutrals:
            if not species.solve_density:
                continue
            
            name = species.name
            
            # Get transport coefficients (or use default for neutrals)
            if name in self.transport.coefficients:
                coeff = self.transport.get(name)
                D = coeff.D
                V = coeff.V
            else:
                # Default neutral diffusion
                D = Function(self.state.V)
                D.x.array[:] = 1.0  # Default D for neutrals
                V = Function(self.state.V)
                V.x.array[:] = 0.0  # No convection for neutrals
            
            # Neutrals always use physics-based diffusion from collision rates
            self._compute_neutral_diffusion(name, D)
            
            # Set source term
            if name in dn_sources:
                self.source.x.array[:] = dn_sources[name]
            else:
                self.source.x.array[:] = 0.0
            
            # Determine BC type
            use_robin = (species.bc_type == "robin")
            bcs = []
            has_dirichlet = (not use_robin and name in self.dirichlet_values)
            n_bc = self.dirichlet_values.get(name, 0.0) if has_dirichlet else 0.0
            
            if has_dirichlet:
                bcs = self.bc_handler.get_dirichlet_bc(n_bc, "both")
            
            # For Robin BC species (H): update decay lengths from physics
            # λ = 2 * D / (vth * (1 - RH))
            if use_robin:
                m_amu = species_mass.get(name, 1.0)
                self._update_robin_decay_lengths(species, D, m_amu)
            
            # Solve density equation
            self.transport_eq.solve(
                D, V, species.n, species.n_prev, species.n_prev2,
                self.source, bcs=bcs, use_robin_bc=use_robin
            )
            
            # Solve energy equation
            # C++ transport.cpp uses gEd, gEv scaling for energy transport:
            # B = Ts - coef1 * gEd * Ds * tstep + coef1 * gEv * Vs * tstep
            if species.solve_energy and name in dE_sources:
                self.source.x.array[:] = dE_sources[name]
                
                # Get neutral energy transport scaling factor (gEdn for neutrals)
                # C++ uses gEdn for neutral energy (boundary conditions), not gEd
                gEdn = self.params.get('gEdn', 5/3)
                
                # Create scaled transport coefficients for energy
                # Neutrals have no convection (V=0), so no gEv needed
                D_energy = Function(self.state.V)
                V_energy = Function(self.state.V)
                D_energy.x.array[:] = gEdn * D.x.array[:]
                V_energy.x.array[:] = 0.0  # Neutrals have no convection
                
                # For Robin BC species (H): use energy-specific decay lengths
                # C++ uses: dE/dr = ±1/2 * gEdn * vth/D * n * 3/2 * T * (1 - REH)
                # This differs from density BC by factor gEdn*(1-REH)/(1-RH)
                if use_robin:
                    m_amu = species_mass.get(name, 1.0)
                    self._update_robin_energy_decay_lengths(species, D, m_amu)
                
                self.transport_eq.solve(
                    D_energy, V_energy, species.E, species.E_prev, species.E_prev2,
                    self.source, bcs=[], use_robin_bc=use_robin
                )
                
                # Apply flux-based energy BC for Dirichlet density species (H2, HeI)
                # - Inward flux: incoming at Ta0, reflected at REH * T_interior
                # - Outward flux: particles leave with their local temperature
                if has_dirichlet:
                    self._apply_energy_bc_with_flux(
                        species.E, species.n, D, n_bc, Ta0, REH
                    )
    
    def _solve_electron_energy(self, dE_sources: dict) -> None:
        """Solve electron energy equation."""
        electrons = self.state.electrons
        
        if not electrons.solve_energy:
            return
        
        # Get electron transport coefficients
        if 'e' in self.transport.coefficients:
            coeff = self.transport.get('e')
            D = coeff.D
            V = coeff.V
        else:
            # Use ion transport as proxy
            D = Function(self.state.V)
            D.x.array[:] = 1.0
            V = Function(self.state.V)
            V.x.array[:] = 0.0
        
        # Set source term
        if 'e' in dE_sources:
            self.source.x.array[:] = dE_sources['e']
        else:
            self.source.x.array[:] = 0.0
        
        # Get energy transport scaling factors (from C++ transport.cpp)
        # C = Ts - coef1 * gEd * Ds * tstep + coef1 * gEv * Vs * tstep
        gEd = self.params.get('gEd', 5/3)
        gEv = self.params.get('gEv', 5/3)
        
        # Create scaled transport coefficients for energy
        D_energy = Function(self.state.V)
        V_energy = Function(self.state.V)
        D_energy.x.array[:] = gEd * D.x.array[:]
        V_energy.x.array[:] = gEv * V.x.array[:]
        
        # Solve with Robin BC
        self.transport_eq.solve(
            D_energy, V_energy, electrons.E, electrons.E_prev, electrons.E_prev2,
            self.source, bcs=[], use_robin_bc=True
        )
    
    def _adapt_timestep(self) -> None:
        """Adapt time step based on solution change."""
        max_change = self.state.max_relative_change()
        
        if max_change > 0:
            # Adjust dt to meet accuracy target
            factor = self.accuracy / max_change
            factor = np.clip(factor, 0.5, 2.0)  # Limit change per step
            
            self.dt = self.dt * factor
            self.dt = np.clip(self.dt, self.dt_min, self.dt_max)
    
    def run_until(self, t_end: float, callback: Callable = None, debug: bool = False) -> None:
        """
        Run simulation until specified time.
        
        Parameters
        ----------
        t_end : float
            End time [s].
        callback : Callable, optional
            Function called after each step: callback(solver, t)
        debug : bool, optional
            Print debug information at each step.
        """
        while self.t < t_end:
            dt = self.step(debug=debug)
            
            if callback is not None:
                callback(self, self.t)


def run_simulation(
    params_or_file,
    output_dir: str = None,
    show_plotter: bool = True
) -> PlasmaState:
    """
    Run a full simulation from input file or parameters dict.
    
    Parameters
    ----------
    params_or_file : str or dict
        Path to JSON input file, or dict of parameters.
    output_dir : str, optional
        Directory for output files.
    show_plotter : bool
        If True, launch interactive Bokeh plotter (default: True).
        
    Returns
    -------
    state : PlasmaState
        Final plasma state.
    """
    from .io.json_input import load_input_file
    from .io.output import write_csv_output
    from .mesh import create_mesh_from_geometry, create_mesh_from_json
    
    # Load parameters
    if isinstance(params_or_file, str):
        params = load_input_file(params_or_file)
    else:
        params = params_or_file
    
    # Create mesh
    if 'grid_file' in params:
        mesh, radial_positions = create_mesh_from_json(params['grid_file'])
    else:
        mesh, radial_positions = create_mesh_from_geometry(
            R=params.get('R', 1.0),
            a=params.get('a', 0.2),
            lHFS=params.get('lHFS', 0.2),
            lLFS=params.get('lLFS', 0.2),
            num_cells=params.get('nmeshp', 100) - 1
        )
    
    # Create plasma state
    state = PlasmaState(mesh)
    
    # Add helium if enabled
    if params.get('bHe', False):
        state.add_helium_species()
    
    # Add molecular hydrogen species (H2, H2+, H3+) if enabled
    if params.get('bH2', False):
        state.add_molecular_hydrogen_species()
    
    # Initialize from parameters
    state.initialize_from_params(params.get('initial_conditions', {}))
    
    # Create boundary condition handler
    bc_handler = BoundaryConditions(mesh, state.V)
    
    # Create transport manager
    transport = TransportManager(state.V)
    transport.initialize_from_params(
        params.get('transport', {}),
        [s.name for s in state.all_species]
    )
    
    # Merge time_step params with main params for limiter access
    solver_params = {**params, **params.get('time_step', {})}
    
    # Create solver with radial positions for limiter losses
    solver = BDF2Solver(state, transport, bc_handler, solver_params, radial_positions)
    
    # Set Dirichlet values for H2, HeI
    if 'nH2_bc' in params:
        solver.set_dirichlet_value('H2', params['nH2_bc'])
    if 'nHeI_bc' in params:
        solver.set_dirichlet_value('HeI', params['nHeI_bc'])
    
    # Output callback - append all outputs to a single file
    output_times = []
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_filename = f"Res_{timestamp}.csv"
    
    # Full path to output file for plotter
    csv_filepath = None
    if output_dir:
        import os
        os.makedirs(output_dir, exist_ok=True)
        csv_filepath = os.path.join(output_dir, output_filename)
    
    # Launch plotter if requested
    plotter_process = None
    if show_plotter and csv_filepath:
        try:
            from .gui.plotter import launch_plotter
            plotter_process = launch_plotter(csv_filepath, show=True)
        except ImportError as e:
            print(f"[Warning] Could not launch plotter: {e}")
            print("[Warning] Install bokeh to enable interactive plotting.")
    
    def output_callback(solver, t):
        output_interval = params.get('output_interval', 1e-4)
        if len(output_times) == 0 or t - output_times[-1] >= output_interval:
            output_times.append(t)
            if output_dir:
                write_csv_output(solver.state, t, radial_positions, output_dir, 
                                filename=output_filename, append=True)
            # Get pecabs from coupled power if available (only for relevant modes)
            pecabs_str = ""
            show_pecabs = (params.get('bgray', False) or params.get('bram', False) or 
                          params.get('bfixpowerfrac', False) or params.get('bnefix', False))
            if show_pecabs and hasattr(solver, 'coupled_power') and solver.coupled_power is not None:
                pecabs_str = f", pecabs = {solver.coupled_power.state.pecabs:.4f}"
            print(f"t = {t:.6e} s, dt = {solver.dt:.6e} s{pecabs_str}")
    
    # Run simulation
    try:
        t_end = params.get('tmainend', 1e-3)
        solver.run_until(t_end, callback=output_callback)
    finally:
        # Stop plotter when simulation ends
        if plotter_process is not None:
            try:
                from .gui.plotter import stop_plotter
                stop_plotter(plotter_process)
            except Exception:
                pass
    
    return state
