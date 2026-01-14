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
from .boundary import BoundaryConditions, DecayLengthBC, compute_neutral_decay_length
from .transport import (TransportCoefficients, TransportManager, 
                        compute_neutral_diffusion_from_state, compute_neutral_diffusion,
                        compute_gyrogeom_diffusion_from_state, DiffusionModel)
from .reactions.collisions import compute_sources_from_state
from .parallel import compute_limiter_losses, compute_bpol_losses, compute_parallel_loss_rates


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
        lambda_n: float = 0.02,
        solver_tolerance: float = 1e-10
    ):
        """
        Initialize transport equation solver.
        
        Parameters
        ----------
        V : FunctionSpace
            Function space for the solution.
        bc_handler : BoundaryConditions
            Boundary condition handler.
        lambda_n : float
            Ion density decay length at boundaries [m] (same for HFS and LFS).
        solver_tolerance : float
            Linear solver tolerance for PETSc KSP.
        """
        self.V = V
        self.mesh = V.mesh
        self.bc_handler = bc_handler
        self.solver_tolerance = solver_tolerance
        
        # Create Robin BC handler (same decay length for HFS and LFS)
        self.robin_bc = DecayLengthBC(lambda_n, lambda_n, bc_handler)
        
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
        
        # ================================================================
        # Persistent Functions for form caching
        # These hold the coefficient values; forms reference them symbolically
        # ================================================================
        self._D_func = fem.Function(V, name="D_cached")
        self._V_func = fem.Function(V, name="V_cached")
        self._n_prev_func = fem.Function(V, name="n_prev_cached")
        self._n_prev2_func = fem.Function(V, name="n_prev2_cached")
        self._source_func = fem.Function(V, name="source_cached")
        self._destruction_rate_func = fem.Function(V, name="kd_cached")
        
        # Cached compiled forms (built on first solve)
        self._a_compiled_robin = None
        self._L_compiled_robin = None
        self._a_compiled_dirichlet = None
        self._L_compiled_dirichlet = None
        
        # Cached KSP solver
        self._ksp = None
        
        # Cached matrix and vector (for reuse)
        self._A_robin = None
        self._A_dirichlet = None
        self._b = None
        self._n_new = None
        
        # Build and cache the forms
        self._build_cached_forms()
    
    def _build_cached_forms(self) -> None:
        """
        Build and compile forms once, using persistent Functions as placeholders.
        
        The form structure never changes - only the Function values change.
        By compiling once and reusing, we save ~80% of solve time.
        """
        r = self.r
        n = self.n_trial
        v = self.v_test
        dt = self.dt
        c1, c2, c3 = self.c1, self.c2, self.c3
        
        # Use persistent Functions as coefficient placeholders
        D = self._D_func
        V_adv = self._V_func
        n_prev = self._n_prev_func
        n_prev2 = self._n_prev2_func
        source = self._source_func
        destruction_rate = self._destruction_rate_func
        
        dx = ufl.dx
        
        # === Build bilinear form (LHS) ===
        # Time derivative: r * (1/dt) * n * v
        # NOTE: Following C++ Tomator formulation where mass term is NOT scaled by c1
        # C++: A = Ts - coef1 * Ds * tstep (mass Ts has no c1, spatial has c1)
        a_time = r * n / dt * v * dx
        
        # Diffusion: c1 * r * D * ∇n · ∇v
        a_diff = c1 * r * D * ufl.inner(ufl.grad(n), ufl.grad(v)) * dx
        
        # Advection: conservative form
        a_adv = -c1 * r * V_adv * n * ufl.grad(v)[0] * dx
        
        # Implicit destruction term: k_d * n on LHS
        a_destruction = c1 * r * destruction_rate * n * v * dx
        
        # Base form without Robin BC
        a_base = a_time + a_diff + a_adv + a_destruction
        
        # === Build linear form (RHS) ===
        L_time = r * (c2 * n_prev - c3 * n_prev2) / dt * v * dx
        L_source = c1 * r * source * v * dx
        L = L_time + L_source
        
        # === Compile forms for Robin BC case ===
        a_boundary = self.robin_bc.get_weak_form_terms(n, v, D, r, c1)
        a_robin = a_base + a_boundary
        
        self._a_compiled_robin = fem.form(a_robin)
        self._L_compiled_robin = fem.form(L)
        
        # === Compile forms for Dirichlet BC case ===
        # Same forms but without Robin boundary terms
        self._a_compiled_dirichlet = fem.form(a_base)
        self._L_compiled_dirichlet = fem.form(L)
        
        # Pre-allocate solution Function
        self._n_new = fem.Function(self.V)
        
        # === Pre-create KSP solver ===
        self._ksp = PETSc.KSP().create(self.mesh.comm)
        self._ksp.setType(PETSc.KSP.Type.PREONLY)
        self._ksp.getPC().setType(PETSc.PC.Type.LU)
        self._ksp.setTolerances(rtol=self.solver_tolerance, atol=self.solver_tolerance, max_it=1000)
        self._ksp.setFromOptions()
    
    def assemble_forms(
        self,
        D: Function,
        V_adv: Function,
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
        V_adv : Function
            Advection velocity [m/s].
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
        # Weak form: (1/dt)*n*v + c1*(diffusion + advection) = (c2/dt)*n_prev - (c3/dt)*n_prev2 + c1*source
        
        # Time derivative: r * (1/dt) * n * v (NOT scaled by c1!)
        a_time = r * n / dt * v * dx
        
        # Diffusion: c1 * r * D * ∇n · ∇v (integration by parts of ∇·(r D ∇n))
        # SCALED by c1 to match C++ Tomator
        a_diff = c1 * r * D * ufl.inner(ufl.grad(n), ufl.grad(v)) * dx
        
        # Advection: CONSERVATIVE weak form of (1/r) ∂/∂r(r V n)
        # Integration by parts: ∫ ∂/∂r(r V n) v dr = -∫ r V n ∂v/∂r dr + [r V n v]_boundary
        # This correctly handles spatially-varying V(r)
        # C++ has: A = Ts - coef1*Ds*tstep + coef1*Vs*tstep
        # So advection has OPPOSITE sign to diffusion on LHS
        # SCALED by c1 to match C++ Tomator
        a_adv = -c1 * r * V_adv * n * ufl.grad(v)[0] * dx
        
        # Combine bilinear form
        a = a_time + a_diff + a_adv
        
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
        V_adv: Function,
        n: Function,
        n_prev: Function,
        n_prev2: Function,
        source: Function,
        bcs: list = None,
        use_robin_bc: bool = True,
        destruction_rate: Function = None
    ) -> None:
        """
        Solve the transport equation for one time step using cached forms.
        
        Parameters
        ----------
        D : Function
            Diffusion coefficient.
        V_adv : Function
            Advection velocity.
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
        # === Copy coefficient values to persistent Functions ===
        # This updates the values used in the cached compiled forms
        self._D_func.x.array[:] = D.x.array[:]
        self._V_func.x.array[:] = V_adv.x.array[:]
        self._n_prev_func.x.array[:] = n_prev.x.array[:]
        self._n_prev2_func.x.array[:] = n_prev2.x.array[:]
        self._source_func.x.array[:] = source.x.array[:]
        
        # Handle optional destruction rate
        if destruction_rate is not None:
            self._destruction_rate_func.x.array[:] = destruction_rate.x.array[:]
        else:
            self._destruction_rate_func.x.array[:] = 0.0
        
        # Select cached forms based on BC type
        if use_robin_bc:
            a_compiled = self._a_compiled_robin
            L_compiled = self._L_compiled_robin
        else:
            a_compiled = self._a_compiled_dirichlet
            L_compiled = self._L_compiled_dirichlet
        
        if bcs is None:
            bcs = []
        
        # === Assemble matrix ===
        A = assemble_matrix(a_compiled, bcs=bcs)
        A.assemble()
        
        # === Assemble RHS vector ===
        b = assemble_vector(L_compiled)
        apply_lifting(b, [a_compiled], [bcs])
        b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
        set_bc(b, bcs)
        
        # === Solve using cached KSP (just update operators) ===
        self._ksp.setOperators(A)
        self._ksp.solve(b, self._n_new.x.petsc_vec)
        self._n_new.x.scatter_forward()
        
        # Copy result to n
        n.x.array[:] = self._n_new.x.array[:]
        
        # Clean up PETSc objects (KSP is reused)
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
            - bc_ion_lambda_n: Ion density decay length [m] (same for HFS/LFS)
            - bc_ion_lambda_E: Ion energy decay length [m]
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
            # Extract DOF coordinates (works for any polynomial degree)
            self.radial_positions = state.V.tabulate_dof_coordinates()[:, 0].copy()
        
        # Time stepping parameters
        self.dt_init = params.get('dtinit', 1e-6)  # Store initial dt
        self.dt = self.dt_init  # Start with dt_init
        self.dt_prev = self.dt
        self.dt_min = params.get('dtmin', 1e-9)
        self.dt_max = params.get('dtmax', 1e-3)
        self.accuracy = params.get('accur', 0.05)
        
        # Step rejection parameters
        self.step_rejection = params.get('step_rejection', True)  # Enable/disable
        self.rejection_margin = params.get('rejection_margin', 1.5)  # Reject if max_change > margin * accur
        self.rejection_safety = params.get('rejection_safety', 0.8)  # Use safety * optimal_dt for retry
        self.max_rejections = params.get('max_rejections', 5)  # Max retries before accepting anyway
        self._rejection_count = 0  # Track rejections for diagnostics
        
        # Ion boundary condition decay length (same for HFS and LFS)
        lambda_n = params.get('bc_ion_lambda_n', 0.02)  # Default 2 cm
        
        # Linear solver tolerance
        solver_tolerance = params.get('solvertolerance', 1e-10)
        
        # Create transport equation solver
        self.transport_eq = TransportEquation(
            state.V, bc_handler, lambda_n, solver_tolerance
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
        
        # Pre-allocate destruction rate Functions for implicit parallel losses
        self._k_n_func = Function(state.V, name="k_n")
        self._k_E_func = Function(state.V, name="k_E")
        
        # Initialize dt to dt_init (will be adapted after each step)
        self.dt = self.dt_init
    
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
        m_amu: float = 1.0,
        is_ion: bool = False
    ) -> None:
        """
        Update Robin BC decay lengths based on local physics (C++ approach).
        
        For neutral H:
            λ = 2 * D / (vth * (1 - R))
        
        For ions (Hi, H2i, H3i, HeII, HeIII):
            λ = sqrt(D * L_conn / (c_s * corrections))
        
        where:
            L_conn = 0.66 * 2π * a_R / n_limiters  [connection length]
            c_s = 9790 * sqrt(Z/μ * (Te + Ti))     [sound speed, m/s]
        
        Parameters
        ----------
        species : Species
            Species to compute decay length for.
        D : Function
            Diffusion coefficient function.
        m_amu : float
            Particle mass in atomic mass units.
        is_ion : bool
            If True, use ion formula with connection length physics.
            If False, use neutral formula with thermal velocity.
        """
        coords = self.state.V.tabulate_dof_coordinates()[:, 0]
        D_arr = D.x.array
        
        # Get temperature at boundaries from energy/density
        E_arr = species.E.x.array
        n_arr = species.n.x.array
        
        sorted_idx = np.argsort(coords)
        
        # HFS boundary (smallest r)
        idx_hfs = sorted_idx[0]
        r_hfs = coords[idx_hfs]
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
        r_lfs = coords[idx_lfs]
        if n_arr[idx_lfs] > 1e-30:
            T_lfs = E_arr[idx_lfs] / (1.5 * n_arr[idx_lfs])
        else:
            # Use interior temperature as fallback
            idx_interior = sorted_idx[-2]
            if n_arr[idx_interior] > 1e-30:
                T_lfs = E_arr[idx_interior] / (1.5 * n_arr[idx_interior])
            else:
                T_lfs = self.params.get('Ta0', 0.026)
        
        if is_ion:
            # Use fixed decay length for density from JSON (C++ lam_dec_length_ions)
            lambda_hfs = self.params.get('bc_ion_lambda_n', 0.02)  # Default 2 cm
            lambda_lfs = self.params.get('bc_ion_lambda_n', 0.02)
        else:
            # Neutral H decay length: thermal velocity formula
            RH = self.params.get('RH', 0.5)
            
            lambda_hfs = compute_neutral_decay_length(
                D=D_arr[idx_hfs],
                T=T_hfs,
                m_amu=m_amu,
                R=RH
            )
            
            lambda_lfs = compute_neutral_decay_length(
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
        m_amu: float = 1.0,
        is_ion: bool = False
    ) -> None:
        """
        Update Robin BC decay lengths for energy equation (C++ approach).
        
        For neutral H:
            The C++ energy BC for H is:
                dE/dr = ±1/2 * gEdn * vth/D * n * 3/2 * T * (1 - RH*REH)
            This uses REH (energy reflection) and includes gEdn factor.
        
        For ions:
            The C++ uses: bEion = bEion / (lambda / sqrt(gEe) * sqrt(gEd))
            So λ_E = λ_n * sqrt(gEe) / sqrt(gEd) = λ_n * sqrt(gEe/gEd)
        
        Parameters
        ----------
        species : Species
            Species to compute decay length for.
        D : Function
            Diffusion coefficient function.
        m_amu : float
            Particle mass in atomic mass units.
        is_ion : bool
            If True, use ion formula with gEe/gEd scaling.
            If False, use neutral formula with gEdn and REH.
        """
        # Get parameters
        RH = self.params.get('RH', 0.5)
        REH = self.params.get('REH', 0.9)
        gEdn = self.params.get('gEdn', 5/3)
        gEd = self.params.get('gEd', 5/3)
        gEe = self.params.get('gEe', 5/3)
        
        coords = self.state.V.tabulate_dof_coordinates()[:, 0]
        D_arr = D.x.array
        
        # Get temperature at boundaries from energy/density
        E_arr = species.E.x.array
        n_arr = species.n.x.array
        
        sorted_idx = np.argsort(coords)
        
        # HFS boundary (smallest r)
        idx_hfs = sorted_idx[0]
        r_hfs = coords[idx_hfs]
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
        r_lfs = coords[idx_lfs]
        if n_arr[idx_lfs] > 1e-30:
            T_lfs = E_arr[idx_lfs] / (1.5 * n_arr[idx_lfs])
        else:
            idx_interior = sorted_idx[-2]
            if n_arr[idx_interior] > 1e-30:
                T_lfs = E_arr[idx_interior] / (1.5 * n_arr[idx_interior])
            else:
                T_lfs = self.params.get('Ta0', 0.026)
        
        if is_ion:
            # Use fixed decay length for energy from JSON (C++ lam_dec_length_energy)
            # In C++, energy uses a DIFFERENT decay length (1 cm) than density (2 cm)
            lambda_E_hfs = self.params.get('bc_ion_lambda_E', 0.01)  # Default 1 cm
            lambda_E_lfs = self.params.get('bc_ion_lambda_E', 0.01)
            
            # Update Robin BC with energy decay lengths directly
            self.transport_eq.robin_bc.update_decay_lengths(lambda_E_hfs, lambda_E_lfs)
        else:
            # Neutral H: use thermal velocity formula with energy corrections
            lambda_n_hfs = compute_neutral_decay_length(
                D=D_arr[idx_hfs],
                T=T_hfs,
                m_amu=m_amu,
                R=RH
            )
            
            lambda_n_lfs = compute_neutral_decay_length(
                D=D_arr[idx_lfs],
                T=T_lfs,
                m_amu=m_amu,
                R=RH
            )
            
            # Update Robin BC with energy-specific decay lengths
            self.transport_eq.robin_bc.update_energy_decay_lengths(
                lambda_n_hfs, lambda_n_lfs, gEdn, REH, RH
            )

    def _compute_ion_gyrogeom_diffusion(self) -> None:
        """
        Compute Gyro-Geometric Diffusion (GGD) coefficient for ion species.
        
        When bDgyrogeom is enabled, computes D using the formula:
            D = Dfsave * (1/3) * nu * mfp * (rho + mfp * Bh/Bt)
        
        This requires collision frequencies (from _nu_collision computed in step()).
        Applies the same D to all ion species: Hi, H2i, H3i, HeII, HeIII.
        
        The GGD model accounts for:
        - Mean free path (mfp) from collisions
        - Gyro radius (rho) from magnetic field (toroidal field Br = Bt * R0 / R)
        - Cross-field transport due to drift (Bh/Bt term)
        """
        # Check if gyrogeom is enabled
        if not self.params.get('bDgyrogeom', False):
            return
        
        # Check if we have collision frequencies
        nu_collision = getattr(self, '_nu_collision', None)
        if nu_collision is None:
            return
        
        # Get magnetic field parameters
        Bt = self.params.get('Bt', 1.0)   # Toroidal field [T]
        Bh = self.params.get('Bh', 0.0)   # Horizontal/poloidal field [T]
        R0 = self.params.get('R', 1.0)    # Major radius [m]
        
        # Dfsave scaling factor
        Dfsave = self.params.get('Dfact', 1.0)
        
        # Compute Br array: toroidal field with 1/R variation
        # C++ uses: Br[id] = Bt * R / aR[id] where R is major radius, aR is radial position
        # Br is the confining toroidal field, NOT a radial component
        R_positions = self.radial_positions  # Radial positions [m]
        B_toroidal = Bt * R0 / R_positions  # [T] - toroidal field with 1/R scaling
        
        # Compute GGD using the dedicated function
        D_gyrogeom = compute_gyrogeom_diffusion_from_state(
            state=self.state,
            nu_collision=nu_collision,
            B_radial=B_toroidal,  # This is actually Br (toroidal) not a radial component
            Bh=Bh,
            Bt=Bt,
            Dfsave=Dfsave,
            D_min=0.01,  # 0.01 m²/s minimum
            D_max=100    # 100 m²/s maximum
        )
        
        # Apply to all ion species that have GYROGEOM model
        ion_species = ['Hi', 'H2i', 'H3i', 'HeII', 'HeIII']
        for name in ion_species:
            if name not in self.transport.coefficients:
                continue
            coeff = self.transport.get(name)
            if coeff.D_model == DiffusionModel.GYROGEOM:
                coeff.D.x.array[:] = D_gyrogeom

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
        coords = self.state.V.tabulate_dof_coordinates()[:, 0]  # DOF coordinates
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
    
    def _compute_max_change(self) -> float:
        """
        Compute maximum relative change across all species (for timestep adaptation).
        
        Ignores vacuum regions (n < nevac) where relative changes can be huge
        but physically irrelevant.
        
        Returns
        -------
        max_change : float
            Maximum relative change in density or energy across all species.
        """
        max_change = 0.0
        nevac = self.params.get('nevac', 1e10)  # Vacuum density threshold [m^-3]
        
        for species in self.state.all_species:
            n_new = species.n.x.array
            n_old = species.n_prev.x.array
            E_new = species.E.x.array
            E_old = species.E_prev.x.array
            
            # Only consider cells above vacuum threshold
            mask = n_new > nevac
            if not np.any(mask):
                continue
            
            with np.errstate(divide='ignore', invalid='ignore'):
                rel_change_n = np.abs(n_new - n_old) / np.maximum(np.abs(n_old), 1e-30)
                rel_change_E = np.abs(E_new - E_old) / np.maximum(np.abs(E_old), 1e-30)
            
            max_change = max(max_change, np.max(rel_change_n[mask]), np.max(rel_change_E[mask]))
        
        # Electrons (energy only) - use electron density for vacuum check
        if self.state.electrons is not None:
            ne = self.state.electrons.n.x.array
            E_new = self.state.electrons.E.x.array
            E_old = self.state.electrons.E_prev.x.array
            
            mask = ne > nevac
            if np.any(mask):
                with np.errstate(divide='ignore', invalid='ignore'):
                    rel_change_E = np.abs(E_new - E_old) / np.maximum(np.abs(E_old), 1e-30)
                max_change = max(max_change, np.max(rel_change_E[mask]))
        
        return max_change

    def _compute_post_solve_dt(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute per-cell optimal dt based on actual changes after solving.
        
        This is for analysis purposes: what would the optimal dt have been
        given the actual changes that occurred in this timestep?
        
        dt = accur × min(n/|Δn|, E/|ΔE|) where Δn = n_new - n_old
        
        Computes separately for:
        - Charged particles: ions (Hi, H2i, H3i, HeII, HeIII) + electrons
        - Neutral particles: H, H2, HeI
        
        Returns
        -------
        dt_charged_total : np.ndarray
            Per-cell optimal dt for all charged particles [s].
        dt_neutral_total : np.ndarray
            Per-cell optimal dt for all neutral particles [s].
        """
        n_points = len(self.radial_positions)
        accur = self.accuracy
        dt = self.dt
        
        # Initialize with large value (no constraint)
        dt_charged = np.full(n_points, 1e30)
        dt_neutral = np.full(n_points, 1e30)
        
        # Ion species
        ion_names = ['Hi', 'H2i', 'H3i', 'HeII', 'HeIII']
        for name in ion_names:
            if name not in self.state.species:
                continue
            species = self.state.species[name]
            
            # Density change
            n_new = species.n.x.array
            n_old = species.n_prev.x.array
            dn = n_new - n_old
            
            # dt = accur * n / |dn/dt| = accur * n * dt / |dn|
            n_safe = np.maximum(np.abs(n_new), 1e-30)
            dn_safe = np.maximum(np.abs(dn), 1e-30)
            dt_n = accur * n_safe * dt / dn_safe
            
            # Only apply where change is significant
            mask = np.abs(dn) > 1e-30
            dt_charged[mask] = np.minimum(dt_charged[mask], dt_n[mask])
            
            # Energy change
            E_new = species.E.x.array
            E_old = species.E_prev.x.array
            dE = E_new - E_old
            
            E_safe = np.maximum(np.abs(E_new), 1e-30)
            dE_safe = np.maximum(np.abs(dE), 1e-30)
            dt_E = accur * E_safe * dt / dE_safe
            
            mask = np.abs(dE) > 1e-30
            dt_charged[mask] = np.minimum(dt_charged[mask], dt_E[mask])
        
        # Electron energy change
        if self.state.electrons is not None:
            E_new = self.state.electrons.E.x.array
            E_old = self.state.electrons.E_prev.x.array
            dE = E_new - E_old
            
            E_safe = np.maximum(np.abs(E_new), 1e-30)
            dE_safe = np.maximum(np.abs(dE), 1e-30)
            dt_E = accur * E_safe * dt / dE_safe
            
            mask = np.abs(dE) > 1e-30
            dt_charged[mask] = np.minimum(dt_charged[mask], dt_E[mask])
        
        # Neutral species
        neutral_names = ['H', 'H2', 'HeI']
        for name in neutral_names:
            if name not in self.state.species:
                continue
            species = self.state.species[name]
            
            # Density change
            n_new = species.n.x.array
            n_old = species.n_prev.x.array
            dn = n_new - n_old
            
            n_safe = np.maximum(np.abs(n_new), 1e-30)
            dn_safe = np.maximum(np.abs(dn), 1e-30)
            dt_n = accur * n_safe * dt / dn_safe
            
            mask = np.abs(dn) > 1e-30
            dt_neutral[mask] = np.minimum(dt_neutral[mask], dt_n[mask])
            
            # Energy change
            E_new = species.E.x.array
            E_old = species.E_prev.x.array
            dE = E_new - E_old
            
            E_safe = np.maximum(np.abs(E_new), 1e-30)
            dE_safe = np.maximum(np.abs(dE), 1e-30)
            dt_E = accur * E_safe * dt / dE_safe
            
            mask = np.abs(dE) > 1e-30
            dt_neutral[mask] = np.minimum(dt_neutral[mask], dt_E[mask])
        
        # Clamp to reasonable range
        dt_charged = np.clip(dt_charged, self.dt_min, self.dt_max)
        dt_neutral = np.clip(dt_neutral, self.dt_min, self.dt_max)
        
        return dt_charged, dt_neutral

    # =========================================================================
    # Helper methods for cleaner code organization
    # =========================================================================

    def _get_species_n(self, name: str) -> Optional[np.ndarray]:
        """Get species density array, or None if species doesn't exist."""
        if name in self.state.species:
            return self.state.species[name].n.x.array.copy()
        return None

    def _get_species_T(self, name: str) -> Optional[np.ndarray]:
        """Get species temperature array, or None if species doesn't exist."""
        if name in self.state.species:
            return self.state.species[name].T.copy()
        return None

    def _get_species_E(self, name: str) -> Optional[np.ndarray]:
        """Get species energy density array [eV·m^-3], or None if species doesn't exist."""
        if name in self.state.species:
            sp = self.state.species[name]
            return 1.5 * sp.T * sp.n.x.array
        return None

    def _merge_sources(self, target: dict, source: dict) -> None:
        """Merge source dict into target dict, adding values for existing keys."""
        for name, arr in source.items():
            if arr is None:
                continue
            if name in target:
                target[name] += arr
            else:
                target[name] = arr.copy()

    # =========================================================================
    # Transport coefficient calculation (equivalent to C++ transpCoef)
    # =========================================================================

    def _compute_transport_coefficients(self) -> None:
        """
        Compute transport coefficients D and V for all species.
        
        Equivalent to C++ transpCoef() function.
        
        Order:
        1. Compute B_toroidal(R) = Bt * R0 / R
        2. Update ion D (fixed, Bohm, or gyrogeom)
        3. Update neutral D (always physics-based from collisions)
        
        Must be called AFTER collision sources are computed (nu_collision available).
        """
        Bt = self.params.get('Bt', 1.0)
        R0 = self.params.get('R', 1.0)
        
        # Compute toroidal field with 1/R scaling: Br = Bt * R0 / R
        # Store for reuse in bpol and other calculations
        self._B_toroidal = Bt * R0 / self.radial_positions
        
        # Update standard transport (fixed or Bohm)
        self.transport.update_all(self.state, B_field=self._B_toroidal)
        
        # Compute gyrogeom diffusion if enabled (uses stored _nu_collision)
        if self.params.get('bDgyrogeom', False):
            self._compute_ion_gyrogeom_diffusion()
        
        # Compute neutral diffusion (always physics-based from collision frequencies)
        for name in ('H', 'H2', 'HeI'):
            if name in self.transport.coefficients:
                coeff = self.transport.get(name)
                self._compute_neutral_diffusion(name, coeff.D)

    # =========================================================================
    # Parallel loss calculations (equivalent to C++ limiters + bpol_function)
    # =========================================================================

    def _add_parallel_losses(self, dn_sources: dict, dE_sources: dict) -> None:
        """
        Add parallel transport losses to source terms.
        
        Includes:
        - Limiter losses (particles hitting poloidal limiters in SOL)
        - Vertical diffusion losses (bpol_function)
        
        Both modify dn_sources and dE_sources in-place.
        """
        # Limiter losses (requires lHFS, lLFS)
        lHFS = self.params.get('lHFS')
        lLFS = self.params.get('lLFS')
        if lHFS is not None and lLFS is not None:
            self._add_limiter_losses(dn_sources, dE_sources)
        
        # Vertical diffusion losses (requires Bt, Bv, b, R)
        if self.params.get('bpol', True) and all(self.params.get(k) is not None for k in ['Bt', 'Bv', 'b', 'R']):
            self._add_bpol_losses(dn_sources, dE_sources)

    def _compute_parallel_loss_rates_k(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute parallel loss rates k_n and k_E for implicit LHS treatment.
        
        Combines limiter losses (SOL region) and bpol losses (vertical diffusion)
        into unified loss rate arrays that are the SAME for all charged species.
        
        Returns
        -------
        k_n : np.ndarray
            Density loss rate 1/τ_n [s^-1].
        k_E : np.ndarray  
            Energy loss rate 1/τ_E [s^-1].
        """
        # Get perpendicular diffusion coefficient (use Hi as representative)
        if 'Hi' in self.transport.coefficients:
            D_perp = self.transport.get('Hi').D.x.array.copy()
        else:
            D_perp = np.ones_like(self.state.electrons.n.x.array) * 1.0
        
        # Determine diffusion model
        if self.params.get('bDgyrogeom', False):
            diffusion_model = 'gyrogeom'
        elif self.params.get('bDbohm', False):
            diffusion_model = 'bohm'
        else:
            diffusion_model = 'fixed'
        
        # Get D_ion for fixed/bohm models
        D_ion = None
        if diffusion_model in ('fixed', 'bohm') and 'Hi' in self.transport.coefficients:
            D_ion = self.transport.get('Hi').D.x.array.copy()
        
        # Get required parameters with defaults
        lHFS = self.params.get('lHFS', self.radial_positions[0])
        lLFS = self.params.get('lLFS', self.radial_positions[-1])
        Br = getattr(self, '_B_toroidal', None)
        Bv = self.params.get('Bv')
        b = self.params.get('b')
        if b is not None:
            b = b * 100.0  # m to cm
        
        k_n, k_E = compute_parallel_loss_rates(
            R_positions=self.radial_positions,
            lHFS=lHFS,
            lLFS=lLFS,
            D_perp=D_perp,
            lambda_n=self.params.get('bc_ion_lambda_n', 0.02),
            lambda_E=self.params.get('bc_ion_lambda_E', 0.01),
            Br=Br,
            Bv=Bv,
            b=b,
            ne=self.state.electrons.n.x.array,
            Te=self.state.electrons.T,
            nHi=self._get_species_n('Hi'),
            THi=self._get_species_T('Hi'),
            nuHi=self._nu_collision.get('Hi') if hasattr(self, '_nu_collision') else None,
            nH2i=self._get_species_n('H2i'),
            TH2i=self._get_species_T('H2i'),
            nuH2i=self._nu_collision.get('H2i') if hasattr(self, '_nu_collision') else None,
            nH3i=self._get_species_n('H3i'),
            TH3i=self._get_species_T('H3i'),
            nuH3i=self._nu_collision.get('H3i') if hasattr(self, '_nu_collision') else None,
            nHeII=self._get_species_n('HeII'),
            THeII=self._get_species_T('HeII'),
            nuHeII=self._nu_collision.get('HeII') if hasattr(self, '_nu_collision') else None,
            nHeIII=self._get_species_n('HeIII'),
            THeIII=self._get_species_T('HeIII'),
            nuHeIII=self._nu_collision.get('HeIII') if hasattr(self, '_nu_collision') else None,
            Dfsave=self.params.get('Dfact', 1.0),
            D_ion=D_ion,
            diffusion_model=diffusion_model,
            gEd=self.params.get('gEd', 1.0)
        )
        
        return k_n, k_E
    def _add_limiter_losses(self, dn_sources: dict, dE_sources: dict) -> None:
        """
        Add limiter parallel losses.
        
        Equivalent to C++ limiters() function.
        Particles in SOL region (r < lHFS or r > lLFS) are lost to limiters.
        Uses diffusion-based formula: τ = λ² / D_⊥
        
        Note: Rate limiting is now applied centrally in step() after all parallel losses computed.
        """
        # Get perpendicular diffusion coefficient (use Hi as representative)
        if 'Hi' in self.transport.coefficients:
            D_perp = self.transport.get('Hi').D.x.array.copy()
        else:
            # Fallback to a reasonable default
            D_perp = np.ones_like(self.state.electrons.n.x.array) * 1.0  # 1 m²/s
        
        dn_lim, dE_lim = compute_limiter_losses(
            R_positions=self.radial_positions,
            lHFS=self.params['lHFS'],
            lLFS=self.params['lLFS'],
            D_perp=D_perp,
            lambda_n=self.params.get('bc_ion_lambda_n', 0.02),  # default 2 cm
            lambda_E=self.params.get('bc_ion_lambda_E', 0.01),  # default 1 cm
            ne=self.state.electrons.n.x.array,
            Te=self.state.electrons.T,
            nHi=self._get_species_n('Hi'),
            THi=self._get_species_T('Hi'),
            nH2i=self._get_species_n('H2i'),
            TH2i=self._get_species_T('H2i'),
            nH3i=self._get_species_n('H3i'),
            TH3i=self._get_species_T('H3i'),
            nHeII=self._get_species_n('HeII'),
            THeII=self._get_species_T('HeII'),
            nHeIII=self._get_species_n('HeIII'),
            THeIII=self._get_species_T('HeIII'),
            Ta0=self.params.get('Ta0', 0.026)
        )
        self._merge_sources(dn_sources, dn_lim)
        self._merge_sources(dE_sources, dE_lim)

    def _add_bpol_losses(self, dn_sources: dict, dE_sources: dict) -> None:
        """
        Add vertical diffusion losses.
        
        Equivalent to C++ bpol_function().
        Vertical losses due to drift in vertical magnetic field component.
        
        For Fixed and Bohm diffusion models, uses the pre-computed radial
        diffusion coefficient. For Gyrogeom, computes vertical-specific
        formula with Bv/Br term.
        """
        ne = self.state.electrons.n.x.array
        Te = self.state.electrons.T
        
        # Determine diffusion model type
        if self.params.get('bDgyrogeom', False):
            diffusion_model = 'gyrogeom'
        elif self.params.get('bDbohm', False):
            diffusion_model = 'bohm'
        else:
            diffusion_model = 'fixed'
        
        # Get pre-computed ion diffusion coefficient (use Hi as representative)
        # For Fixed/Bohm, this is used directly. For Gyrogeom, it's ignored.
        D_ion = None
        if diffusion_model in ('fixed', 'bohm'):
            if 'Hi' in self.transport.coefficients:
                D_ion = self.transport.get('Hi').D.x.array.copy()
        
        dn_bpol, dE_bpol = compute_bpol_losses(
            Br=self._B_toroidal,
            Bv=self.params['Bv'],
            b=self.params['b'] * 100.0,  # m to cm for CGS formulas
            ne=ne,
            Te=Te,
            nHi=self._get_species_n('Hi'),
            THi=self._get_species_T('Hi'),
            EHi=self._get_species_E('Hi'),
            nuHi=self._nu_collision.get('Hi', np.ones_like(ne) * 5e3),
            nH2i=self._get_species_n('H2i'),
            TH2i=self._get_species_T('H2i'),
            EH2i=self._get_species_E('H2i'),
            nuH2i=self._nu_collision.get('H2i'),
            nH3i=self._get_species_n('H3i'),
            TH3i=self._get_species_T('H3i'),
            EH3i=self._get_species_E('H3i'),
            nuH3i=self._nu_collision.get('H3i'),
            nHeII=self._get_species_n('HeII'),
            THeII=self._get_species_T('HeII'),
            EHeII=self._get_species_E('HeII'),
            nuHeII=self._nu_collision.get('HeII'),
            nHeIII=self._get_species_n('HeIII'),
            THeIII=self._get_species_T('HeIII'),
            EHeIII=self._get_species_E('HeIII'),
            nuHeIII=self._nu_collision.get('HeIII'),
            Dfsave=self.params.get('Dfact', 1.0),
            D_ion=D_ion,
            diffusion_model=diffusion_model,
            gEd=self.params.get('gEd', 1.0),
            Ta0=self.params.get('Ta0', 0.026)
        )
        self._merge_sources(dn_sources, dn_bpol)
        self._merge_sources(dE_sources, dE_bpol)

    def _solve_neutrals_transport_only(self, dn_sources: dict, dE_sources: dict) -> None:
        """
        Solve transport equations for neutral species (transport only).
        
        For operator splitting: no reaction sources, only diffusion/advection.
        """
        Ta0 = self.params.get('Ta0', 0.026)
        REH = self.params.get('REH', 0.9)
        species_mass = {'H': 1.0, 'H2': 2.0, 'HeI': 4.0}
        
        for species in self.state.neutrals:
            if not species.solve_density:
                continue
            
            name = species.name
            
            if name in self.transport.coefficients:
                coeff = self.transport.get(name)
                D = coeff.D
                V = coeff.V
            else:
                D = Function(self.state.V)
                D.x.array[:] = 1.0
                V = Function(self.state.V)
                V.x.array[:] = 0.0
            
            self._compute_neutral_diffusion(name, D)
            
            # Parallel losses only (no reaction sources for neutrals in transport step)
            if name in dn_sources:
                self.source.x.array[:] = dn_sources[name]
            else:
                self.source.x.array[:] = 0.0
            
            use_robin = (species.bc_type == "robin")
            bcs = []
            has_dirichlet = (not use_robin and name in self.dirichlet_values)
            n_bc = self.dirichlet_values.get(name, 0.0) if has_dirichlet else 0.0
            
            if has_dirichlet:
                bcs = self.bc_handler.get_dirichlet_bc(n_bc, "both")
            
            if use_robin:
                m_amu = species_mass.get(name, 1.0)
                self._update_robin_decay_lengths(species, D, m_amu, is_ion=False)
            
            self.transport_eq.solve(
                D, V, species.n, species.n_prev, species.n_prev2,
                self.source, bcs=bcs, use_robin_bc=use_robin
            )
            
            # Solve energy equation (always, sources are handled by reaction step)
            if species.solve_energy:
                if name in dE_sources:
                    self.source.x.array[:] = dE_sources[name]
                else:
                    self.source.x.array[:] = 0.0
                
                gEdn = self.params.get('gEdn', 5/3)
                D_energy = Function(self.state.V)
                V_energy = Function(self.state.V)
                D_energy.x.array[:] = gEdn * D.x.array[:]
                V_energy.x.array[:] = 0.0
                
                if use_robin:
                    m_amu = species_mass.get(name, 1.0)
                    self._update_robin_energy_decay_lengths(species, D, m_amu, is_ion=False)
                
                self.transport_eq.solve(
                    D_energy, V_energy, species.E, species.E_prev, species.E_prev2,
                    self.source, bcs=[], use_robin_bc=use_robin
                )

    # =========================================================================
    # RF power deposition
    # =========================================================================

    def _add_rf_power(self, dE_sources: dict) -> None:
        """
        Add RF power deposition to electron energy source.
        
        Equivalent to C++ coupled power calculation.
        """
        if not hasattr(self, 'coupled_power') or self.coupled_power is None:
            return
        
        result = self.coupled_power.compute_power(
            R=self.radial_positions,
            ne=self.state.electrons.n.x.array,
            Te=self.state.electrons.T,
            dt=self.dt,
            t=self.t,
            nue=self._nu_collision.get('e')
        )
        
        # Add RF power to electron energy
        if 'e' in dE_sources:
            dE_sources['e'] += result.PRFe
        else:
            dE_sources['e'] = result.PRFe.copy()
        
        # Store absorbed power fraction for diagnostics
        self._pecabs = result.pecabs

    # =========================================================================
    # Post-solve corrections
    # =========================================================================

    def _recompute_electron_density_preserve_temperature(self) -> None:
        """
        Recompute ne from quasi-neutrality while preserving Te.
        
        After limiting ion solutions, ne must be recomputed but Te should
        stay constant (not rescaled with ne).
        """
        ne_before = self.state.electrons.n.x.array.copy()
        self.state.compute_electron_density()
        ne_after = self.state.electrons.n.x.array
        
        # Scale energy to preserve temperature: E_new = E_old * (n_new / n_old)
        scale = np.where(ne_before > 1e-30, ne_after / ne_before, 1.0)
        self.state.electrons.E.x.array[:] *= scale

    def _recompute_energy_from_temperature(self) -> None:
        """
        Recompute E = 1.5 * n * T for all species to ensure consistency.
        
        Called after temperature clamps to guarantee E and T are consistent.
        This is the final step of post-solve corrections.
        """
        ENERGY_FACTOR = 1.5
        
        # All heavy species (ions and neutrals)
        for species in self.state.all_species:
            n_arr = species.n.x.array
            T_arr = species.T  # T property computes E / (1.5 * n)
            species.E.x.array[:] = ENERGY_FACTOR * n_arr * T_arr
        
        # Electrons
        if self.state.electrons is not None:
            n_arr = self.state.electrons.n.x.array
            T_arr = self.state.electrons.T
            self.state.electrons.E.x.array[:] = ENERGY_FACTOR * n_arr * T_arr

    # =========================================================================
    # Timestep calculation methods
    # =========================================================================

    def _compute_dt_collisions(
        self, 
        dn_sources: Dict[str, np.ndarray], 
        dE_sources: Dict[str, np.ndarray]
    ) -> np.ndarray:
        """
        Compute per-cell timestep limit from collision + RF sources.
        
        For each grid point, computes the dt that would make the maximum
        relative change equal to accur, considering all species (n and E).
        
        dt[i] = accur × min(n/|dn|, E/|dE|) across all species at point i
        
        Grid points where a species has n < nevac are ignored for that species
        (large relative changes in vacuum regions don't impact physics).
        
        Parameters
        ----------
        dn_sources : dict
            Density source terms [m^-3/s] for each species.
        dE_sources : dict
            Energy source terms [eV·m^-3/s] for each species.
            
        Returns
        -------
        dt_collision : np.ndarray
            Per-cell timestep limit [s].
        """
        n_points = len(self.radial_positions)
        
        # Initialize with large value (no constraint)
        dt_limit = np.full(n_points, 1e30)
        
        accur = self.accuracy
        nevac = self.params.get('nevac', 1e10)  # Vacuum density threshold
        
        # Check all species densities
        for name, dn in dn_sources.items():
            if dn is None:
                continue
            
            # Get current density
            if name == 'e':
                n = self.state.electrons.n.x.array
            elif name in self.state.species:
                n = self.state.species[name].n.x.array
            else:
                continue
            
            # dt = accur * n / |dn| where dn != 0
            n_safe = np.maximum(np.abs(n), 1e-30)
            dn_safe = np.maximum(np.abs(dn), 1e-30)
            dt_n = accur * n_safe / dn_safe
            
            # Only apply where dn is significant AND density is above vacuum
            # (ignore vacuum regions where relative changes are large but irrelevant)
            mask = (np.abs(dn) > 1e-30) & (n > nevac)
            dt_limit[mask] = np.minimum(dt_limit[mask], dt_n[mask])
        
        # Check all species energies
        for name, dE in dE_sources.items():
            if dE is None:
                continue
            
            # Get current energy and density
            if name == 'e':
                E = self.state.electrons.E.x.array
                n = self.state.electrons.n.x.array
            elif name in self.state.species:
                E = self.state.species[name].E.x.array
                n = self.state.species[name].n.x.array
            else:
                continue
            
            # dt = accur * E / |dE| where dE != 0
            E_safe = np.maximum(np.abs(E), 1e-30)
            dE_safe = np.maximum(np.abs(dE), 1e-30)
            dt_E = accur * E_safe / dE_safe
            
            # Only apply where dE is significant AND density is above vacuum
            mask = (np.abs(dE) > 1e-30) & (n > nevac)
            dt_limit[mask] = np.minimum(dt_limit[mask], dt_E[mask])
        
        # Clamp to reasonable range
        dt_limit = np.clip(dt_limit, self.dt_min, self.dt_max)
        
        return dt_limit

    def _compute_dt_diffusion(self) -> Dict[str, np.ndarray]:
        """
        Compute per-cell diffusion stability timestep limit.
        
        Uses the stability criterion: dt ≤ safety × Δx² / (2D)
        with safety factor 0.5.
        
        Computes separately for:
        - ions: collective ion diffusion coefficient
        - H: neutral hydrogen
        - H2: molecular hydrogen
        - HeI: neutral helium
        
        Returns
        -------
        dt_diffusion : dict
            Per-cell timestep limits for each species group.
            Keys: 'ions', 'H', 'H2', 'HeI'
        """
        safety = 0.5
        
        # Compute per-cell Δx (distance between adjacent nodes)
        dx = np.diff(self.radial_positions)  # length n_points - 1
        
        # Extend to n_points by duplicating last value
        dx_cell = np.zeros(len(self.radial_positions))
        dx_cell[:-1] = dx
        dx_cell[-1] = dx[-1]
        
        dt_diffusion = {}
        
        # Ion diffusion (use Hi as representative)
        if 'Hi' in self.transport.coefficients:
            D_ion = self.transport.get('Hi').D.x.array
            D_safe = np.maximum(D_ion, 1e-10)
            dt_diffusion['ions'] = safety * dx_cell**2 / (2.0 * D_safe)
        
        # Neutral species
        for name in ('H', 'H2', 'HeI'):
            if name in self.transport.coefficients:
                D_neut = self.transport.get(name).D.x.array
                D_safe = np.maximum(D_neut, 1e-10)
                dt_diffusion[name] = safety * dx_cell**2 / (2.0 * D_safe)
        
        # Clamp all to reasonable range
        for key in dt_diffusion:
            dt_diffusion[key] = np.clip(dt_diffusion[key], self.dt_min, self.dt_max)
        
        return dt_diffusion

    def _compute_optimal_dt(
        self,
        dt_collision: np.ndarray,
        dt_diffusion: Dict[str, np.ndarray]
    ) -> float:
        """
        Compute optimal timestep combining collision and diffusion limits.
        
        Uses harmonic mean: dt = (1/dt_col + 1/dt_diff)^{-1}
        This ensures the smaller limit dominates.
        
        Takes the global minimum across all grid points.
        
        Parameters
        ----------
        dt_collision : np.ndarray
            Per-cell collision timestep limit [s].
        dt_diffusion : dict
            Per-cell diffusion timestep limits for each species group.
            
        Returns
        -------
        dt_optimal : float
            Optimal global timestep [s].
        """
        n_points = len(self.radial_positions)
        
        # Start with collision limit
        inv_dt = 1.0 / np.maximum(dt_collision, 1e-30)
        
        # Add ion diffusion limit only (ignore neutral diffusion for optimal dt)
        # Neutrals have much larger D (~1000x ions) so their dt_diff would be overly restrictive
        if 'ions' in dt_diffusion:
            inv_dt += 1.0 / np.maximum(dt_diffusion['ions'], 1e-30)
        
        # Harmonic mean: dt = 1 / (sum of 1/dt_i)
        dt_per_cell = 1.0 / inv_dt
        
        # Global minimum
        dt_optimal = np.min(dt_per_cell)
        
        # Clamp to bounds (but allow dt < dt_min on first step to grow from dt_init)
        if self.step_count > 0:
            dt_optimal = np.clip(dt_optimal, self.dt_min, self.dt_max)
        else:
            dt_optimal = min(dt_optimal, self.dt_max)
        
        return dt_optimal

    def _save_state_for_rejection(self) -> dict:
        """
        Save current state so we can restore if step is rejected.
        
        Returns
        -------
        saved_state : dict
            Dictionary containing copies of all state arrays.
        """
        saved = {}
        
        # Save all species (ions and neutrals)
        for species in self.state.all_species:
            name = species.name
            saved[f'{name}_n'] = species.n.x.array.copy()
            saved[f'{name}_E'] = species.E.x.array.copy()
        
        # Save electrons
        saved['e_n'] = self.state.electrons.n.x.array.copy()
        saved['e_E'] = self.state.electrons.E.x.array.copy()
        
        return saved
    
    def _restore_state_from_saved(self, saved: dict) -> None:
        """
        Restore state from saved copy after step rejection.
        
        Parameters
        ----------
        saved : dict
            Dictionary containing saved state arrays from _save_state_for_rejection.
        """
        # Restore all species (ions and neutrals)
        for species in self.state.all_species:
            name = species.name
            species.n.x.array[:] = saved[f'{name}_n']
            species.E.x.array[:] = saved[f'{name}_E']
        
        # Restore electrons
        self.state.electrons.n.x.array[:] = saved['e_n']
        self.state.electrons.E.x.array[:] = saved['e_E']

    def _solve_with_current_dt(
        self,
        k_n: np.ndarray,
        k_E: np.ndarray,
        dn_sources: Dict[str, np.ndarray],
        dE_sources: Dict[str, np.ndarray],
        debug: bool = False
    ) -> None:
        """
        Solve transport equations with current dt.
        
        This is the core solve step that can be retried with different dt.
        Sources and transport coefficients are reused (depend only on previous state).
        
        Parameters
        ----------
        k_n : np.ndarray
            Density loss rate for implicit parallel losses.
        k_E : np.ndarray
            Energy loss rate for implicit parallel losses.
        dn_sources : dict
            Density source terms.
        dE_sources : dict
            Energy source terms.
        debug : bool
            Print debug information.
        """
        # Debug helper
        def _dbg(msg, arrays=None):
            if not debug:
                return
            print(f"  [DBG] {msg}")
            if arrays:
                for name, arr in arrays.items():
                    if arr is None:
                        continue
                    arr_np = arr if isinstance(arr, np.ndarray) else arr.x.array
                    print(f"        {name}: min={arr_np.min():.3e}, max={arr_np.max():.3e}")
        
        # Update BDF2 coefficients with current timestep
        self.transport_eq.update_bdf2_coefficients(self.dt, self.dt_prev, self.step_count)
        
        # Solve ions
        self._solve_ions_with_implicit_k(k_n, k_E, dn_sources, dE_sources)
        
        if debug:
            _dbg("After ion transport:", {
                'nHi': self._get_species_n('Hi'),
                'nHeII': self._get_species_n('HeII'),
            })
        
        # Solve neutrals
        self._solve_neutrals_transport_only(dn_sources, dE_sources)
        
        if debug:
            _dbg("After neutral transport:", {
                'nH': self._get_species_n('H'),
            })
        
        # Solve electron energy
        self._solve_electron_energy_with_implicit_k(k_E, dE_sources)
        
        if debug:
            _dbg("After electron energy transport:", {'Te': self.state.electrons.T})
        
        # Post-solve corrections
        self._apply_density_floors()
        self._recompute_electron_density_preserve_temperature()
        self._apply_temperature_clamps()
        
        if debug:
            _dbg("After post-solve corrections:", {'ne': self.state.electrons.n, 'Te': self.state.electrons.T})

    # =========================================================================
    # Main time stepping method
    # =========================================================================

    def step(self, debug: bool = False, profile: bool = False) -> float:
        """
        Perform one time step with adaptive timestep control and optional step rejection.
        
        Timestep strategy:
        - dt grows or shrinks so that max_change stays near accur
        - After each step: dt_next = dt * accur / max_change (clamped to dt_max)
        - Optional step rejection: if max_change > rejection_margin * accur,
          restore state and retry with smaller dt
        
        Calculation order:
        1. Compute collision sources and collision frequencies (nu)
        2. Compute RF power coupling
        3. Compute collision timestep limit (dt_collision) - for analysis
        4. Compute parallel loss rates k_n, k_E for implicit treatment
        5. Compute transport coefficients D, V
        6. Compute diffusion timestep limit (dt_diffusion) - for analysis
        7. Solve transport with step rejection loop:
           - Save state, solve, check max_change
           - If max_change > rejection_margin * accur: restore and retry with smaller dt
           - Sources/transport coefficients reused (depend on previous state only)
        8. Adapt dt for next step based on final max_change
        
        Parameters
        ----------
        debug : bool
            Print debug information.
        profile : bool
            Print timing breakdown for each step section.
        
        Returns
        -------
        dt : float
            The time step that was taken.
        """
        import time
        timings = {}
        
        def tic(name):
            if profile:
                timings[name] = time.perf_counter()
        
        def toc(name):
            if profile:
                timings[name] = time.perf_counter() - timings[name]
        
        # Debug helper
        def _dbg(msg, arrays=None):
            if not debug:
                return
            print(f"  [DBG] {msg}")
            if arrays:
                for name, arr in arrays.items():
                    if arr is None:
                        continue
                    arr_np = arr if isinstance(arr, np.ndarray) else arr.x.array
                    print(f"        {name}: min={arr_np.min():.3e}, max={arr_np.max():.3e}")
        
        # Timestep was already adapted at end of previous step
        # (or set to dt_init in __init__)
        
        if debug:
            print(f"\n=== STEP {self.step_count}, t={self.t:.3e}, dt={self.dt:.3e} ===")
            _dbg("Initial state:", {
                'ne': self.state.electrons.n,
                'nHi': self._get_species_n('Hi'),
                'nHeII': self._get_species_n('HeII'),
            })
        
        # ================================================================
        # STEP 1: Compute collision sources and frequencies (nu)
        # ================================================================
        tic('collisions')
        dn_sources, dE_sources, nu_collision = compute_sources_from_state(
            self.state, self.params
        )
        self._nu_collision = nu_collision
        
        # Apply energy factor: E = 3/2 * n * T
        ENERGY_FACTOR = 1.5
        for species_name in dE_sources:
            dE_sources[species_name] *= ENERGY_FACTOR
        toc('collisions')
        
        # ================================================================
        # STEP 2: Compute RF power coupling
        # ================================================================
        tic('rf_power')
        dE_RF = {}
        self._add_rf_power(dE_RF)
        
        # Merge RF into collision sources for timestep calculation
        self._merge_sources(dE_sources, dE_RF)
        toc('rf_power')
        
        # ================================================================
        # STEP 3: Compute collision timestep limit (per-cell)
        # dt_col = accur × min(n/|dn|, E/|dE|) across all species
        # ================================================================
        tic('dt_collision')
        dt_collision = self._compute_dt_collisions(dn_sources, dE_sources)
        self._dt_collision = dt_collision  # Store for output
        toc('dt_collision')
        
        # ================================================================
        # STEP 4: Compute parallel loss rates k_n, k_E for implicit treatment
        # Now that nu is available, compute the shared loss rates
        # ================================================================
        tic('parallel_loss_rates')
        # First need transport coefficients for D_perp
        # Compute transport coefficients early (needed for parallel loss rates)
        self._compute_transport_coefficients()
        
        # Compute k_n and k_E (same for all charged species)
        k_n, k_E = self._compute_parallel_loss_rates_k()
        
        # Store in Functions for passing to transport solver
        self._k_n_func.x.array[:] = k_n
        self._k_E_func.x.array[:] = k_E
        toc('parallel_loss_rates')
        
        # ================================================================
        # STEP 5: Transport coefficients already computed above
        # ================================================================
        
        # ================================================================
        # STEP 6: Compute diffusion timestep limit (per-cell) - for analysis
        # dt_diff = 0.5 × Δx² / (2D) for each species group
        # ================================================================
        tic('dt_diffusion')
        dt_diffusion = self._compute_dt_diffusion()
        self._dt_diffusion = dt_diffusion  # Store for output
        toc('dt_diffusion')
        
        # ================================================================
        # STEP 7: Compute optimal timestep from collision/diffusion (for analysis)
        # The actual dt used is dt_optimal_prev from previous step
        # ================================================================
        tic('optimal_dt')
        dt_optimal_sources = self._compute_optimal_dt(dt_collision, dt_diffusion)
        self._dt_optimal = dt_optimal_sources  # Store for diagnostics (source-based)
        
        if debug:
            print(f"  [TIMESTEP] Using dt={self.dt:.3e} (dt_optimal_sources={dt_optimal_sources:.3e})")
        toc('optimal_dt')
        
        # ================================================================
        # STEP 8-10: Solve with step rejection
        # Save state, solve, check max_change, reject if too large
        # Sources and transport coefficients are reused (depend on previous state)
        # ================================================================
        tic('solve_with_rejection')
        
        # Save state for potential rejection
        saved_state = self._save_state_for_rejection()
        
        rejection_count = 0
        while True:
            # Solve transport + post-corrections with current dt
            self._solve_with_current_dt(k_n, k_E, dn_sources, dE_sources, debug=debug)
            
            # Check max_change
            max_change = self._compute_max_change()
            
            if debug:
                print(f"  [TIMESTEP] max_change={max_change:.3e}, target={self.accuracy:.3e}, margin={self.rejection_margin:.1f}")
            
            # Check if step should be rejected
            reject_threshold = self.accuracy * self.rejection_margin
            should_reject = (
                self.step_rejection and 
                max_change > reject_threshold and 
                rejection_count < self.max_rejections
            )
            
            if should_reject:
                rejection_count += 1
                self._rejection_count += 1
                
                # Compute optimal dt with safety factor
                dt_optimal = self.dt * (self.accuracy / max_change) * self.rejection_safety
                dt_optimal = max(dt_optimal, self.dt_min * 0.1)  # Don't go too small
                
                if debug:
                    print(f"  [REJECT #{rejection_count}] max_change={max_change:.3e} > {reject_threshold:.3e}")
                    print(f"  [REJECT] Retrying with dt={dt_optimal:.3e} (was {self.dt:.3e})")
                
                # Restore state and retry with smaller dt
                self._restore_state_from_saved(saved_state)
                self.dt = dt_optimal
            else:
                # Step accepted
                if debug and rejection_count > 0:
                    print(f"  [ACCEPTED] after {rejection_count} rejection(s)")
                break
        
        toc('solve_with_rejection')
        
        dt_taken = self.dt

        # ================================================================
        # Compute post-solve dt diagnostics (for analysis)
        # ================================================================
        tic('post_solve_dt')
        dt_charged_total, dt_neutral_total = self._compute_post_solve_dt()
        self._dt_charged_total = dt_charged_total
        self._dt_neutral_total = dt_neutral_total
        toc('post_solve_dt')

        # ================================================================
        # STEP 11: Adapt timestep for next iteration
        # Simple scheme: dt grows/shrinks so max_change stays near accur
        # max_change was already computed in the rejection loop
        # ================================================================
        tic('timestep_adapt')
        
        # Adapt dt for next step: dt_next = dt * accur / max_change
        if max_change > 1e-30:
            dt_next = self.dt * (self.accuracy / max_change)
        else:
            # No change - allow dt to grow to max
            dt_next = self.dt_max
        
        # Clamp to dt_max always, but only enforce dt_min once we've grown above it
        # This allows starting with dt_init < dt_min and growing from there
        dt_next = min(dt_next, self.dt_max)
        if self.dt >= self.dt_min:
            # Once above dt_min, don't go below it
            dt_next = max(dt_next, self.dt_min)
        self.dt = dt_next
        
        if debug:
            print(f"  [TIMESTEP] dt_next={self.dt:.3e}")
        
        toc('timestep_adapt')
        
        # Store previous solutions for BDF2
        self.state.store_all_previous()
        
        # Update time and counter
        self.t += dt_taken
        self.step_count += 1
        self.dt_prev = dt_taken
        
        # Print timing breakdown if profiling enabled
        if profile:
            total = sum(timings.values())
            print(f"  [PROFILE] Step {self.step_count-1}: total={total*1000:.2f}ms")
            for name, t in sorted(timings.items(), key=lambda x: -x[1]):
                pct = 100 * t / total if total > 0 else 0
                print(f"    {name:20s}: {t*1000:7.2f}ms ({pct:5.1f}%)")
        
        return dt_taken
    
    def _solve_ions_with_implicit_k(
        self, 
        k_n: np.ndarray, 
        k_E: np.ndarray,
        dn_sources: Dict[str, np.ndarray],
        dE_sources: Dict[str, np.ndarray]
    ) -> None:
        """
        Solve transport equations for ions with implicit parallel loss rates
        and explicit reaction sources.
        
        The parallel losses (k_n * n and k_E * E) are added to the LHS of the
        transport equation as destruction terms, making them unconditionally stable.
        Reaction sources are added as explicit RHS terms.
        
        Parameters
        ----------
        k_n : np.ndarray
            Density loss rate 1/τ_n [s^-1].
        k_E : np.ndarray
            Energy loss rate 1/τ_E [s^-1].
        dn_sources : dict
            Density source terms [m^-3/s] for each species.
        dE_sources : dict
            Energy source terms [eV·m^-3/s] for each species.
        """
        ion_mass = {'Hi': 1.0, 'H2i': 2.0, 'H3i': 3.0, 'HeII': 4.0, 'HeIII': 4.0}
        
        # Create destruction rate Functions
        k_n_func = Function(self.state.V)
        k_E_func = Function(self.state.V)
        k_n_func.x.array[:] = k_n
        k_E_func.x.array[:] = k_E
        
        for species in self.state.ions:
            if not species.solve_density:
                continue
            
            name = species.name
            
            coeff = self.transport.get(name)
            D = coeff.D
            V = coeff.V
            
            # Reaction sources as explicit RHS term
            if name in dn_sources:
                self.source.x.array[:] = dn_sources[name]
            else:
                self.source.x.array[:] = 0.0
            
            use_robin = (species.bc_type == "robin")
            bcs = []
            
            if not use_robin and name in self.dirichlet_values:
                bcs = self.bc_handler.get_dirichlet_bc(
                    self.dirichlet_values[name], "both"
                )
            
            if use_robin:
                m_amu = ion_mass.get(name, 1.0)
                self._update_robin_decay_lengths(species, D, m_amu, is_ion=True)
            
            # Solve density with implicit k_n
            self.transport_eq.solve(
                D, V, species.n, species.n_prev, species.n_prev2,
                self.source, bcs=bcs, use_robin_bc=use_robin,
                destruction_rate=k_n_func  # Implicit parallel losses
            )
            
            self._enforce_outflow_boundary(species.n)
            
            # Solve energy equation with implicit k_E and reaction sources
            if species.solve_energy:
                if name in dE_sources:
                    self.source.x.array[:] = dE_sources[name]
                else:
                    self.source.x.array[:] = 0.0
                
                gEd = self.params.get('gEd', 5/3)
                gEv = self.params.get('gEv', 5/3)
                
                D_energy = Function(self.state.V)
                V_energy = Function(self.state.V)
                D_energy.x.array[:] = gEd * D.x.array[:]
                V_energy.x.array[:] = gEv * V.x.array[:]
                
                if use_robin:
                    m_amu = ion_mass.get(name, 1.0)
                    self._update_robin_energy_decay_lengths(species, D, m_amu, is_ion=True)
                
                self.transport_eq.solve(
                    D_energy, V_energy, species.E, species.E_prev, species.E_prev2,
                    self.source, bcs=[], use_robin_bc=use_robin,
                    destruction_rate=k_E_func  # Implicit parallel losses
                )
                self._enforce_outflow_boundary(species.E)
    
    def _solve_electron_energy_with_implicit_k(
        self, 
        k_E: np.ndarray, 
        dE_sources: Dict[str, np.ndarray]
    ) -> None:
        """
        Solve electron energy equation with implicit parallel losses and sources.
        
        Sources include both RF power and collision sources (already merged).
        
        Parameters
        ----------
        k_E : np.ndarray
            Energy loss rate 1/τ_E [s^-1].
        dE_sources : dict
            Energy source terms [eV·m^-3/s] including RF power for electrons.
        """
        # Create destruction rate Function
        k_E_func = Function(self.state.V)
        k_E_func.x.array[:] = k_E
        
        electrons = self.state.electrons
        
        # Electron energy sources (collisions + RF, already merged)
        if 'e' in dE_sources:
            self.source.x.array[:] = dE_sources['e']
        else:
            self.source.x.array[:] = 0.0
        
        # Use ion diffusion coefficient scaled by gEe for electrons
        gEe = self.params.get('gEe', 5/3)
        
        if 'Hi' in self.transport.coefficients:
            D_ion = self.transport.get('Hi').D.x.array
        else:
            D_ion = np.ones_like(electrons.n.x.array) * 1.0
        
        D_energy = Function(self.state.V)
        V_energy = Function(self.state.V)
        D_energy.x.array[:] = gEe * D_ion
        V_energy.x.array[:] = 0.0  # No advection for electron energy
        
        # Update Robin BC for electron energy
        lambda_E = self.params.get('bc_ion_lambda_E', 0.01)
        self.transport_eq.robin_bc.update_decay_lengths(lambda_E, lambda_E)
        
        self.transport_eq.solve(
            D_energy, V_energy, electrons.E, electrons.E_prev, electrons.E_prev2,
            self.source, bcs=[], use_robin_bc=True,
            destruction_rate=k_E_func  # Implicit parallel losses
        )
    
    def _apply_density_floors(self) -> None:
        """
        Apply density floor clamping to all ion species while preserving temperature.
        
        This matches C++ Tomator1D timeStep.cpp:
        - nevac = 1.0 cm⁻³ = 1e6 m⁻³ is the floor for all ion densities
        - When density is below floor, scale E to preserve T: E *= nevac/n
        
        Neutrals are NOT floored (they can be very dilute in vacuum regions).
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
                # Scale E to preserve T: E_new = E_old * (nevac / n_old)
                # This matches C++: Er.EHi[im] = Er.EHi[im] * nevac / nr.nHi[im]
                E_arr[below_floor] *= nevac_si / n_arr[below_floor]
                
                # Clamp density to floor
                n_arr[below_floor] = nevac_si
    
    def _apply_temperature_clamps(self) -> None:
        """
        Clamp energy density to enforce temperature limits.
        
        This is CRITICAL for numerical stability when densities are low.
        The temperature T = E / (1.5 * n) can blow up if n is small.
        We enforce E <= 1.5 * n * T_max for all species.
        
        Temperature limits (matching C++ Tomator1D):
        - Heavy particles (ions, neutrals): T_max = 1000 eV
        - Electrons: T_max = 1000 eV
        """
        T_MIN = 0.026  # Room temperature
        T_MAX_HEAVY = 1000.0  # Max for ions/neutrals
        T_MAX_ELECTRON = 1000.0  # Max for electrons (matching C++)
        
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
        coords = self.state.V.tabulate_dof_coordinates()[:, 0]  # DOF coordinates
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
    
    def _adapt_timestep(self, max_change: float = None) -> None:
        """
        Adapt time step based on solution change.
        
        Uses C++ Tomator approach: allow larger dt increases when solution
        is smooth (factor up to maxtstepincrement), but allow unlimited
        decreases when accuracy requires it.
        
        Parameters
        ----------
        max_change : float, optional
            Maximum relative change to use for adaptation. If not provided,
            computes from current state. Pass the RAW change before limiting
            for proper timestep adaptation.
        """
        if max_change is None:
            max_change = self.state.max_relative_change()
        
        # Get maximum time step increment factor from params (default 2.0, C++ uses up to 10)
        max_increment = self.params.get('maxtstepincrement', 2.0)
        
        if max_change > 0:
            # Adjust dt to meet accuracy target
            factor = self.accuracy / max_change
            # Only limit increases (up to maxtstepincrement), allow unlimited decreases
            factor = min(factor, max_increment)
            
            self.dt = self.dt * factor
            # While dt < dt_min, let it grow naturally with maxtstepincrement
            # Once dt >= dt_min, enforce the dt_min floor
            if self.dt >= self.dt_min:
                self.dt = np.clip(self.dt, self.dt_min, self.dt_max)
            else:
                # Below dt_min: only enforce dt_max, let dt grow toward dt_min
                self.dt = min(self.dt, self.dt_max)
    
    def run_until(self, t_end: float, callback: Callable = None, debug: bool = False, profile: bool = False) -> None:
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
        profile : bool, optional
            Print timing breakdown for each step.
        """
        while self.t < t_end:
            dt = self.step(debug=debug, profile=profile)
            
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
    
    # FEM polynomial degree (affects stencil width)
    # degree=1: linear elements (3-point stencil)
    # degree=2: quadratic elements (5-point stencil)
    # degree=3: cubic elements (7-point stencil)
    fem_degree = params.get('fem_degree', 1)
    
    # Target DOFs from nmeshp, adjust mesh cells to maintain ~constant DOFs
    # For degree=d with N cells: DOFs = d*N + 1
    # To get target_dofs with degree d: N = (target_dofs - 1) / d
    target_dofs = params.get('nmeshp', 161)
    num_cells = max(1, (target_dofs - 1) // fem_degree)
    
    # Create mesh
    if 'grid_file' in params:
        mesh, radial_positions = create_mesh_from_json(params['grid_file'])
    else:
        mesh, radial_positions = create_mesh_from_geometry(
            R=params.get('R', 1.0),
            a=params.get('a', 0.2),
            lHFS=params.get('lHFS', 0.2),
            lLFS=params.get('lLFS', 0.2),
            num_cells=num_cells
        )
    
    # Create plasma state
    state = PlasmaState(mesh, degree=fem_degree)
    actual_dofs = state.V.dofmap.index_map.size_global
    print(f"[FEM] degree={fem_degree}, cells={num_cells}, DOFs={actual_dofs} (target={target_dofs})")
    
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
    
    # Create solver - radial_positions will be computed from DOF coordinates
    # (don't pass mesh vertices as they don't match DOFs for degree > 1)
    solver = BDF2Solver(state, transport, bc_handler, solver_params)
    
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
                # Use DOF coordinates from solver (matches function values)
                dof_coords = solver.radial_positions
                
                # Gather transport data - ions and neutrals separately
                transport_data = {}
                
                # Ion diffusion (D_ion, V_ion) - use Hi as representative
                if 'Hi' in solver.transport.coefficients:
                    coeff = solver.transport.get('Hi')
                    transport_data['D'] = coeff.D.x.array.copy()
                    transport_data['V'] = coeff.V.x.array.copy()
                elif 'e' in solver.transport.coefficients:
                    coeff = solver.transport.get('e')
                    transport_data['D'] = coeff.D.x.array.copy()
                    transport_data['V'] = coeff.V.x.array.copy()
                
                # Neutral diffusion coefficients (D_H, D_H2, D_HeI)
                for neutral_name in ['H', 'H2', 'HeI']:
                    if neutral_name in solver.transport.coefficients:
                        coeff = solver.transport.get(neutral_name)
                        transport_data[f'D_{neutral_name}'] = coeff.D.x.array.copy()
                
                if not transport_data:
                    transport_data = None
                
                # Gather coupled power data
                power_data = None
                if hasattr(solver, 'coupled_power') and solver.coupled_power is not None:
                    # Recompute power profile for output
                    from .reactions.collisions import compute_sources_from_state
                    _, _, nu_collision = compute_sources_from_state(solver.state, params)
                    power_result = solver.coupled_power.compute_power(
                        R=dof_coords,
                        ne=solver.state.electrons.n.x.array.copy(),
                        Te=solver.state.electrons.T.copy(),
                        dt=solver.dt,
                        t=t,
                        nue=nu_collision.get('e', None)
                    )
                    power_data = power_result.PRFe
                
                # Gather timestep limit profiles
                dt_data = None
                if hasattr(solver, '_dt_collision') and solver._dt_collision is not None:
                    dt_data = {}
                    dt_data['dt_collision'] = solver._dt_collision.copy()
                    if hasattr(solver, '_dt_diffusion') and solver._dt_diffusion is not None:
                        # Ion diffusion limit
                        if 'ions' in solver._dt_diffusion:
                            dt_data['dt_ion_diff'] = solver._dt_diffusion['ions'].copy()
                        # Neutral diffusion limit (use H as representative, or take min)
                        neutral_keys = [k for k in solver._dt_diffusion if k in ('H', 'H2', 'HeI')]
                        if neutral_keys:
                            # Take minimum across all neutral species
                            dt_neutral = solver._dt_diffusion[neutral_keys[0]].copy()
                            for k in neutral_keys[1:]:
                                dt_neutral = np.minimum(dt_neutral, solver._dt_diffusion[k])
                            dt_data['dt_neutral_diff'] = dt_neutral
                    # Post-solve dt diagnostics (based on actual changes)
                    if hasattr(solver, '_dt_charged_total') and solver._dt_charged_total is not None:
                        dt_data['dt_charged_total'] = solver._dt_charged_total.copy()
                    if hasattr(solver, '_dt_neutral_total') and solver._dt_neutral_total is not None:
                        dt_data['dt_neutral_total'] = solver._dt_neutral_total.copy()
                
                write_csv_output(solver.state, t, dof_coords, output_dir, 
                                filename=output_filename, append=True,
                                transport_data=transport_data, power_data=power_data,
                                dt_data=dt_data)
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
        output_params = params.get('output_parameters', {})
        profile = params.get('profile', output_params.get('profile', False))
        solver.run_until(t_end, callback=output_callback, profile=profile)
    finally:
        # Stop plotter when simulation ends
        if plotter_process is not None:
            try:
                from .gui.plotter import stop_plotter
                stop_plotter(plotter_process)
            except Exception:
                pass
    
    return state
