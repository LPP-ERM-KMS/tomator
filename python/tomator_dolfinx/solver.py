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
from .boundary import BoundaryConditions, DecayLengthBC, compute_neutral_decay_length, compute_ion_decay_length
from .transport import (TransportCoefficients, TransportManager, 
                        compute_neutral_diffusion_from_state, compute_neutral_diffusion,
                        compute_gyrogeom_diffusion_from_state, DiffusionModel)
from .reactions.collisions import compute_sources_from_state
from .reactions.implicit_reactions import solve_reactions_vectorized
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
        # Time derivative: r * (1/dt) * n * v (NOT scaled by c1!)
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
            # Extract DOF coordinates (works for any polynomial degree)
            self.radial_positions = state.V.tabulate_dof_coordinates()[:, 0].copy()
        
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
            # Check if fixed BC is enabled for ions (C++ fixBCs)
            fixBC = self.params.get('fixBC', False)
            
            if fixBC:
                # Use fixed decay length from JSON (converted to meters)
                lambda_hfs = self.params.get('fixBC_value', 0.02)  # Default 2 cm
                lambda_lfs = self.params.get('fixBC_value', 0.02)
            else:
                # Ion decay length: connection length physics
                nlimiters = self.params.get('nlimiters', 1)
                nevac = self.params.get('nevac', 1.0) * 1e6  # Convert cm^-3 to m^-3
                
                # Get electron temperature and density from electrons species
                Te_arr = self.state.electrons.T  # T property returns array [eV]
                ne_arr = self.state.electrons.n.x.array  # n is a fem.Function
                Te_hfs = max(Te_arr[idx_hfs], 0.01)
                Te_lfs = max(Te_arr[idx_lfs], 0.01)
                ne_hfs = max(ne_arr[idx_hfs], 1e10)
                ne_lfs = max(ne_arr[idx_lfs], 1e10)
                
                # Charge state (Z=1 for H+, H2+, H3+, HeII; Z=2 for HeIII)
                name = species.name
                Z = 2.0 if name == 'HeIII' else 1.0
                
                lambda_hfs = compute_ion_decay_length(
                    D=D_arr[idx_hfs],
                    Te=Te_hfs,
                    Ti=T_hfs,
                    a_R=r_hfs,
                    nlimiters=nlimiters,
                    Z=Z,
                    mu=m_amu,
                    ne=ne_hfs,
                    ni=max(n_arr[idx_hfs], 1e10),
                    nevac=nevac
                )
                
                lambda_lfs = compute_ion_decay_length(
                    D=D_arr[idx_lfs],
                    Te=Te_lfs,
                    Ti=T_lfs,
                    a_R=r_lfs,
                    nlimiters=nlimiters,
                    Z=Z,
                    mu=m_amu,
                    ne=ne_lfs,
                    ni=max(n_arr[idx_lfs], 1e10),
                    nevac=nevac
                )
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
            # Check if fixed BC is enabled for ions (C++ fixBCs)
            fixBC = self.params.get('fixBC', False)
            
            if fixBC:
                # Use fixed decay length from JSON (converted to meters)
                # C++ uses lam_dec_length_energy = 1.0 cm for energy, but we use same value
                lambda_n_hfs = self.params.get('fixBC_value', 0.02)  # Default 2 cm
                lambda_n_lfs = self.params.get('fixBC_value', 0.02)
            else:
                # Ion energy decay length: uses connection length formula with gEe/gEd scaling
                # C++: bEion = bEion / (lambda / sqrt(gEe) * sqrt(gEd))
                # So λ_E = λ_n * sqrt(gEe/gEd)
                nlimiters = self.params.get('nlimiters', 1)
                nevac = self.params.get('nevac', 1.0) * 1e6
                
                # Get electron temperature and density from electrons species
                Te_arr = self.state.electrons.T  # T property returns array [eV]
                ne_arr = self.state.electrons.n.x.array  # n is a fem.Function
                Te_hfs = max(Te_arr[idx_hfs], 0.01)
                Te_lfs = max(Te_arr[idx_lfs], 0.01)
                ne_hfs = max(ne_arr[idx_hfs], 1e10)
                ne_lfs = max(ne_arr[idx_lfs], 1e10)
                
                name = species.name
                Z = 2.0 if name == 'HeIII' else 1.0
                
                # Compute density decay length first
                lambda_n_hfs = compute_ion_decay_length(
                    D=D_arr[idx_hfs],
                    Te=Te_hfs,
                    Ti=T_hfs,
                    a_R=r_hfs,
                    nlimiters=nlimiters,
                    Z=Z,
                    mu=m_amu,
                    ne=ne_hfs,
                    ni=max(n_arr[idx_hfs], 1e10),
                    nevac=nevac
                )
                
                lambda_n_lfs = compute_ion_decay_length(
                    D=D_arr[idx_lfs],
                    Te=Te_lfs,
                    Ti=T_lfs,
                    a_R=r_lfs,
                    nlimiters=nlimiters,
                    Z=Z,
                    mu=m_amu,
                    ne=ne_lfs,
                    ni=max(n_arr[idx_lfs], 1e10),
                    nevac=nevac
                )
            
            # Apply gEe/gEd scaling for energy: λ_E = λ_n * sqrt(gEe/gEd)
            scale = np.sqrt(gEe / gEd) if gEd > 0 else 1.0
            lambda_E_hfs = lambda_n_hfs * scale
            lambda_E_lfs = lambda_n_lfs * scale
            
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
    
    def _limit_solution_change(self) -> Tuple[float, float]:
        """
        Limit solution change to prevent extreme jumps.
        
        Only limits changes that are extremely large (> LIMIT_THRESHOLD).
        Normal changes within accuracy target are allowed - timestep adaptation
        handles those. This prevents numerical blowup while not slowing convergence.
        
        Returns
        -------
        max_change_before : float
            The maximum relative change BEFORE limiting (for timestep adaptation).
        max_change_after : float
            The maximum relative change after limiting.
        """
        # Only limit changes larger than this threshold (prevents blowup without slowing convergence)
        LIMIT_THRESHOLD = 0.5  # 50% change limit
        
        max_change_before = 0.0
        max_change_after = 0.0
        
        # Apply per-point limiting to all species
        for species in self.state.all_species:
            n_new = species.n.x.array
            n_old = species.n_prev.x.array
            E_new = species.E.x.array
            E_old = species.E_prev.x.array
            
            # Compute per-point relative change for density
            with np.errstate(divide='ignore', invalid='ignore'):
                rel_change_n = np.abs(n_new - n_old) / np.maximum(np.abs(n_old), 1e-30)
            
            # Track max change BEFORE limiting
            max_change_before = max(max_change_before, np.max(rel_change_n))
            
            # Find points exceeding LIMIT_THRESHOLD (not accuracy!)
            exceed_n = rel_change_n > LIMIT_THRESHOLD
            if np.any(exceed_n):
                # Compute per-point limiting factor (avoid div by zero)
                safe_rel_change_n = np.maximum(rel_change_n, 1e-30)
                factor_n = LIMIT_THRESHOLD / safe_rel_change_n
                # Apply: limited = old + factor * (new - old)
                n_new[exceed_n] = n_old[exceed_n] + factor_n[exceed_n] * (n_new[exceed_n] - n_old[exceed_n])
                max_change_after = max(max_change_after, LIMIT_THRESHOLD)
            else:
                max_change_after = max(max_change_after, np.max(rel_change_n) if len(rel_change_n) > 0 else 0.0)
            
            # Compute per-point relative change for energy
            with np.errstate(divide='ignore', invalid='ignore'):
                rel_change_E = np.abs(E_new - E_old) / np.maximum(np.abs(E_old), 1e-30)
            
            # Track max change BEFORE limiting
            max_change_before = max(max_change_before, np.max(rel_change_E))
            
            # Find points exceeding LIMIT_THRESHOLD
            exceed_E = rel_change_E > LIMIT_THRESHOLD
            if np.any(exceed_E):
                # Compute per-point limiting factor (avoid div by zero)
                safe_rel_change_E = np.maximum(rel_change_E, 1e-30)
                factor_E = LIMIT_THRESHOLD / safe_rel_change_E
                # Apply: limited = old + factor * (new - old)
                E_new[exceed_E] = E_old[exceed_E] + factor_E[exceed_E] * (E_new[exceed_E] - E_old[exceed_E])
                max_change_after = max(max_change_after, LIMIT_THRESHOLD)
            else:
                max_change_after = max(max_change_after, np.max(rel_change_E) if len(rel_change_E) > 0 else 0.0)
        
        # Apply to electrons - only limit ENERGY, not density
        # Electron density is determined by quasi-neutrality, not directly limited
        if self.state.electrons is not None:
            E_new = self.state.electrons.E.x.array
            E_old = self.state.electrons.E_prev.x.array
            
            # Energy limiting only
            with np.errstate(divide='ignore', invalid='ignore'):
                rel_change_E = np.abs(E_new - E_old) / np.maximum(np.abs(E_old), 1e-30)
            
            max_change_before = max(max_change_before, np.max(rel_change_E))
            
            exceed_E = rel_change_E > LIMIT_THRESHOLD
            if np.any(exceed_E):
                safe_rel_change_E = np.maximum(rel_change_E, 1e-30)
                factor_E = LIMIT_THRESHOLD / safe_rel_change_E
                E_new[exceed_E] = E_old[exceed_E] + factor_E[exceed_E] * (E_new[exceed_E] - E_old[exceed_E])
                max_change_after = max(max_change_after, LIMIT_THRESHOLD)
            else:
                max_change_after = max(max_change_after, np.max(rel_change_E) if len(rel_change_E) > 0 else 0.0)
        
        return max_change_before, max_change_after

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
        if all(self.params.get(k) is not None for k in ['Bt', 'Bv', 'b', 'R']):
            self._add_bpol_losses(dn_sources, dE_sources)

    def _add_limiter_losses(self, dn_sources: dict, dE_sources: dict) -> None:
        """
        Add limiter parallel losses.
        
        Equivalent to C++ limiters() function.
        Particles in SOL region (r < lHFS or r > lLFS) are lost to limiters.
        """
        dn_lim, dE_lim = compute_limiter_losses(
            R_positions=self.radial_positions,
            lHFS=self.params['lHFS'],
            lLFS=self.params['lLFS'],
            nlimiters=self.params.get('nlimiters', 6),
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

    def _compute_parallel_losses(self, dn_sources: dict, dE_sources: dict) -> None:
        """
        Compute parallel transport losses without adding to existing sources.
        
        Same as _add_parallel_losses but for operator splitting where we need
        parallel losses separate from reaction sources.
        
        Parameters
        ----------
        dn_sources : dict
            Empty dict, will be populated with parallel loss sources.
        dE_sources : dict
            Empty dict, will be populated with parallel energy loss sources.
        """
        # Limiter losses (requires lHFS, lLFS)
        lHFS = self.params.get('lHFS')
        lLFS = self.params.get('lLFS')
        if lHFS is not None and lLFS is not None:
            dn_lim, dE_lim = compute_limiter_losses(
                R_positions=self.radial_positions,
                lHFS=self.params['lHFS'],
                lLFS=self.params['lLFS'],
                nlimiters=self.params.get('nlimiters', 6),
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
        
        # Vertical diffusion losses (requires Bt, Bv, b, R)
        if all(self.params.get(k) is not None for k in ['Bt', 'Bv', 'b', 'R']):
            # Determine diffusion model type
            if self.params.get('bDgyrogeom', False):
                diffusion_model = 'gyrogeom'
            elif self.params.get('bDbohm', False):
                diffusion_model = 'bohm'
            else:
                diffusion_model = 'fixed'
            
            D_ion = None
            if diffusion_model in ('fixed', 'bohm'):
                if 'Hi' in self.transport.coefficients:
                    D_ion = self.transport.get('Hi').D.x.array.copy()
            
            ne = self.state.electrons.n.x.array
            Te = self.state.electrons.T
            
            dn_bpol, dE_bpol = compute_bpol_losses(
                Br=self._B_toroidal,
                Bv=self.params['Bv'],
                b=self.params['b'] * 100.0,
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
        
        # Apply energy factor: E = 3/2 * n * T
        ENERGY_FACTOR = 1.5
        for species_name in dE_sources:
            dE_sources[species_name] *= ENERGY_FACTOR

    def _solve_ions_transport_only(self, dn_sources: dict, dE_sources: dict) -> None:
        """
        Solve transport equations for ion species (transport + parallel losses only).
        
        For operator splitting: no reaction sources, only diffusion/advection/parallel losses.
        """
        ion_mass = {'Hi': 1.0, 'H2i': 2.0, 'H3i': 3.0, 'HeII': 4.0, 'HeIII': 4.0}
        
        for species in self.state.ions:
            if not species.solve_density:
                continue
            
            name = species.name
            
            coeff = self.transport.get(name)
            D = coeff.D
            V = coeff.V
            
            # Parallel losses only (no reaction sources)
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
            
            self.transport_eq.solve(
                D, V, species.n, species.n_prev, species.n_prev2,
                self.source, bcs=bcs, use_robin_bc=use_robin
            )
            
            self._enforce_outflow_boundary(species.n)
            
            # Solve energy equation
            if species.solve_energy and name in dE_sources:
                self.source.x.array[:] = dE_sources[name]
                
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
                    self.source, bcs=[], use_robin_bc=use_robin
                )
                self._enforce_outflow_boundary(species.E)

    def _solve_neutrals_transport_only(self, dn_sources: dict, dE_sources: dict) -> None:
        """
        Solve transport equations for neutral species (transport + parallel losses only).
        
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
            
            if species.solve_energy and name in dE_sources:
                self.source.x.array[:] = dE_sources[name]
                
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

    # =========================================================================
    # Main time stepping method
    # =========================================================================

    def step(self, debug: bool = False, profile: bool = False) -> float:
        """
        Perform one time step.
        
        The calculation order matches C++ Tomator1D:
        1. Compute collision sources and collision frequencies (nu)
        2. Compute transport coefficients (D, V) - these use nu
        3. Add parallel losses (limiters, bpol) - these use nu and D
        4. Add RF power deposition
        5. Solve transport equations for all species
        6. Apply floors/clamps and update electron density
        7. Solve electron energy
        8. Limit solution change and adapt timestep
        
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
        
        if debug:
            print(f"\n=== STEP {self.step_count}, t={self.t:.3e}, dt={self.dt:.3e} ===")
            _dbg("Initial state:", {
                'ne': self.state.electrons.n,
                'nHi': self._get_species_n('Hi'),
                'nHeII': self._get_species_n('HeII'),
            })
        
        # Check if using operator splitting for reactions
        use_operator_splitting = self.params.get('operator_splitting', True)
        
        # ================================================================
        # STEP 1: Update BDF2 time discretization coefficients
        # ================================================================
        tic('bdf2_coef')
        self.transport_eq.update_bdf2_coefficients(self.dt, self.dt_prev, self.step_count)
        toc('bdf2_coef')
        
        # ================================================================
        # STEP 2: Compute collision frequencies (nu)
        # nu is needed by: transport coefficients, bpol losses, RF power
        # For operator splitting, we only need nu here (sources computed in reaction step)
        # ================================================================
        tic('collisions')
        dn_sources, dE_sources, nu_collision = compute_sources_from_state(self.state, self.params)
        toc('collisions')
        
        # Store nu for transport coefficient and loss calculations
        self._nu_collision = nu_collision
        
        # Apply energy factor: E = 3/2 * n * T, so dE = 3/2 * d(nT)
        ENERGY_FACTOR = 1.5
        for species_name in dE_sources:
            dE_sources[species_name] *= ENERGY_FACTOR
        
        if debug:
            _dbg("Collision sources (dn/dt):", {
                k: v for k, v in dn_sources.items() if v is not None
            })
        
        # ================================================================
        # STEP 3: Compute transport coefficients D and V
        # Equivalent to C++ transpCoef() function
        # Must be done BEFORE bpol because bpol uses D (via Dfsave)
        # ================================================================
        tic('transport_coef')
        self._compute_transport_coefficients()
        toc('transport_coef')
        
        # ================================================================
        # STEP 4: Compute parallel loss sources (limiter + bpol)
        # ================================================================
        tic('parallel_losses')
        dn_parallel = {}
        dE_parallel = {}
        self._compute_parallel_losses(dn_parallel, dE_parallel)
        toc('parallel_losses')
        
        # ================================================================
        # STEP 5: Add RF power deposition to electron energy
        # ================================================================
        tic('rf_power')
        self._add_rf_power(dE_sources)
        toc('rf_power')
        
        # ================================================================
        # STEP 6: TRANSPORT STEP (PDE with parallel losses only)
        # For operator splitting: transport has no reaction sources
        # Parallel losses are included as explicit sources in transport
        # ================================================================
        if use_operator_splitting:
            # Transport with parallel losses only (no reaction sources)
            tic('solve_ions')
            self._solve_ions_transport_only(dn_parallel, dE_parallel)
            toc('solve_ions')
            
            if debug:
                _dbg("After ion transport:", {
                    'nHi': self._get_species_n('Hi'),
                    'nHeII': self._get_species_n('HeII'),
                })
            
            tic('solve_neutrals')
            self._solve_neutrals_transport_only(dn_parallel, dE_parallel)
            toc('solve_neutrals')
            
            if debug:
                _dbg("After neutral transport:", {
                    'nH': self._get_species_n('H'),
                })
            
            # ================================================================
            # STEP 7: REACTION STEP (implicit ODE at each mesh point)
            # ================================================================
            tic('reactions_implicit')
            solve_reactions_vectorized(self.state, self.dt, self.params)
            toc('reactions_implicit')
            
            if debug:
                _dbg("After implicit reactions:", {
                    'nHi': self._get_species_n('Hi'),
                    'nH': self._get_species_n('H'),
                })
        else:
            # Original explicit treatment (all sources combined)
            self._add_parallel_losses(dn_sources, dE_sources)
            
            tic('solve_ions')
            self._solve_ions(dn_sources, dE_sources)
            toc('solve_ions')
            
            if debug:
                _dbg("After ion solve:", {
                    'nHi': self._get_species_n('Hi'),
                    'nHeII': self._get_species_n('HeII'),
                })
            
            tic('solve_neutrals')
            self._solve_neutrals(dn_sources, dE_sources)
            toc('solve_neutrals')
            
            if debug:
                _dbg("After neutral solve:", {
                    'nH': self._get_species_n('H'),
                })
        
        # ================================================================
        # STEP 7: Apply density floors and temperature clamps
        # CRITICAL for stability (matching C++)
        # ================================================================
        tic('floors_clamps')
        self._apply_density_floors()
        self._apply_temperature_clamps()
        toc('floors_clamps')
        
        # Update electron density from quasi-neutrality
        tic('electron_density')
        self.state.compute_electron_density()
        toc('electron_density')
        
        if debug:
            _dbg("After QN ne update:", {'ne': self.state.electrons.n})
        
        # ================================================================
        # STEP 8: Solve electron energy equation
        # ================================================================
        tic('solve_electron_E')
        self._solve_electron_energy(dE_sources)
        toc('solve_electron_E')
        
        if debug:
            _dbg("After energy solve:", {'Te': self.state.electrons.T})
        
        # ================================================================
        # STEP 9: Limit solution change and adapt timestep
        # ================================================================
        tic('limit_adapt')
        max_change_raw, max_change_limited = self._limit_solution_change()
        
        # Recompute ne after ion limiting, preserve Te
        self._recompute_electron_density_preserve_temperature()

        if debug and max_change_limited >= self.accuracy:
            print(f"  [LIMITED] max_change {max_change_raw:.3e} -> {max_change_limited:.3e}")
        
        # Adapt timestep based on raw change
        dt_taken = self.dt
        self._adapt_timestep(max_change_raw)
        
        # Store previous solutions for BDF2
        self.state.store_all_previous()
        toc('limit_adapt')
        
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
    
    def _apply_density_floors(self) -> None:
        """
        Apply density floor clamping to all ion species.
        
        This matches C++ Tomator1D timeStep.cpp:
        - nevac = 1.0 cm⁻³ = 1e6 m⁻³ is the floor for all ion densities
        - When density is below floor, set temperature to vacuum temperature Ta0
        """
        # nevac from C++: 1.0 cm⁻³ = 1e6 m⁻³ in SI units
        nevac_si = self.params.get('nevac', 1.0) * 1e6  # Convert cm⁻³ to m⁻³
        
        # Vacuum temperature (Ta0) for low-density regions
        Ta0 = self.params.get('Ta0', 0.026)  # [eV]
        
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
                # Set energy such that T = Ta0 at low-density points
                # E = 1.5 * n * T, so E = 1.5 * nevac * Ta0
                E_arr[below_floor] = 1.5 * nevac_si * Ta0
                
                # Clamp density to floor
                n_arr[below_floor] = nevac_si
        
        # Also apply floor to neutral species (optional, for stability)
        neutral_names = ['H', 'H2', 'HeI']
        for name in neutral_names:
            if name not in self.state.species:
                continue
            species = self.state.species[name]
            n_arr = species.n.x.array
            E_arr = species.E.x.array
            
            # Use smaller floor for neutrals since they can be very dilute
            # but still set temperature to Ta0 where density is very low
            neutral_floor = 1e-10
            below_floor = n_arr < neutral_floor
            if np.any(below_floor):
                E_arr[below_floor] = 1.5 * neutral_floor * Ta0
                n_arr[below_floor] = neutral_floor
    
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
            # Ions use connection length formula: λ = sqrt(D * L_conn / c_s)
            if use_robin:
                m_amu = ion_mass.get(name, 1.0)
                self._update_robin_decay_lengths(species, D, m_amu, is_ion=True)
            
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
                # Ions use connection length formula with gEe/gEd scaling
                if use_robin:
                    m_amu = ion_mass.get(name, 1.0)
                    self._update_robin_energy_decay_lengths(species, D, m_amu, is_ion=True)
                
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
                V.x.array[:] = 0.0  # No advection for neutrals
            
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
            # Neutral H uses thermal velocity formula: λ = 2D / (vth * (1 - RH))
            if use_robin:
                m_amu = species_mass.get(name, 1.0)
                self._update_robin_decay_lengths(species, D, m_amu, is_ion=False)
            
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
                # Neutrals have no advection (V=0), so no gEv needed
                D_energy = Function(self.state.V)
                V_energy = Function(self.state.V)
                D_energy.x.array[:] = gEdn * D.x.array[:]
                V_energy.x.array[:] = 0.0  # Neutrals have no advection
                
                # For Robin BC species (H): use energy-specific decay lengths
                # C++ uses: dE/dr = ±1/2 * gEdn * vth/D * n * 3/2 * T * (1 - REH)
                # This differs from density BC by factor gEdn*(1-REH)/(1-RH)
                if use_robin:
                    m_amu = species_mass.get(name, 1.0)
                    self._update_robin_energy_decay_lengths(species, D, m_amu, is_ion=False)
                
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
            self.dt = np.clip(self.dt, self.dt_min, self.dt_max)
    
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
                
                write_csv_output(solver.state, t, dof_coords, output_dir, 
                                filename=output_filename, append=True,
                                transport_data=transport_data, power_data=power_data)
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
