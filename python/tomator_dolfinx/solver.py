"""
Main solver for 1D cylindrical plasma transport using dolfinx.

Implements:
- Weak form assembly for reaction-diffusion-advection equations
- BDF2 adaptive time stepping
- Species transport with Robin/Dirichlet boundary conditions
"""

from typing import Optional, Dict, Tuple, Callable
import numpy as np
import ufl
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import fem, default_scalar_type
from dolfinx.fem import Function, FunctionSpace, Constant
from dolfinx.fem.petsc import assemble_matrix, assemble_vector, apply_lifting, set_bc

from .species import Species, PlasmaState
from .boundary import BoundaryConditions, DecayLengthBC
from .transport import TransportCoefficients, TransportManager
from .reactions.sources import compute_sources_from_state


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
        
        # Solve with KSP
        ksp = PETSc.KSP().create(self.mesh.comm)
        ksp.setOperators(A)
        ksp.setType(PETSc.KSP.Type.GMRES)
        ksp.getPC().setType(PETSc.PC.Type.ILU)
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
        params: dict
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
        """
        self.state = state
        self.transport = transport
        self.bc_handler = bc_handler
        self.params = params
        
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
        
        # 2. Compute collision source terms
        dn_sources, dE_sources = compute_sources_from_state(self.state, self.params)
        
        if debug:
            _dbg("Source terms (dn/dt):", {
                k: v for k, v in dn_sources.items() if v is not None
            })
        
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
        for species in self.state.ions:
            if not species.solve_density:
                continue
            
            name = species.name
            
            # Get transport coefficients
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
            
            # Solve density equation
            self.transport_eq.solve(
                D, V, species.n, species.n_prev, species.n_prev2,
                self.source, bcs=bcs, use_robin_bc=use_robin
            )
            
            # CRITICAL: Enforce outflow-only at boundaries for ions
            # This prevents artificial ion influx that causes numerical blowup
            self._enforce_outflow_boundary(species.n)
            
            # Solve energy equation (similar structure)
            if species.solve_energy and name in dE_sources:
                self.source.x.array[:] = dE_sources[name]
                self.transport_eq.solve(
                    D, V, species.E, species.E_prev, species.E_prev2,
                    self.source, bcs=[], use_robin_bc=use_robin
                )
                # Also enforce outflow for energy
                self._enforce_outflow_boundary(species.E)
    
    def _solve_neutrals(self, dn_sources: dict, dE_sources: dict) -> None:
        """Solve transport equations for neutral species."""
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
            
            # Set source term
            if name in dn_sources:
                self.source.x.array[:] = dn_sources[name]
            else:
                self.source.x.array[:] = 0.0
            
            # Determine BC type
            use_robin = (species.bc_type == "robin")
            bcs = []
            
            if not use_robin and name in self.dirichlet_values:
                bcs = self.bc_handler.get_dirichlet_bc(
                    self.dirichlet_values[name], "both"
                )
            
            # Solve density equation
            self.transport_eq.solve(
                D, V, species.n, species.n_prev, species.n_prev2,
                self.source, bcs=bcs, use_robin_bc=use_robin
            )
            
            # Solve energy equation
            if species.solve_energy and name in dE_sources:
                self.source.x.array[:] = dE_sources[name]
                self.transport_eq.solve(
                    D, V, species.E, species.E_prev, species.E_prev2,
                    self.source, bcs=[], use_robin_bc=use_robin
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
        
        # Solve with Robin BC
        self.transport_eq.solve(
            D, V, electrons.E, electrons.E_prev, electrons.E_prev2,
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
    output_dir: str = None
) -> PlasmaState:
    """
    Run a full simulation from input file or parameters dict.
    
    Parameters
    ----------
    params_or_file : str or dict
        Path to JSON input file, or dict of parameters.
    output_dir : str, optional
        Directory for output files.
        
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
    
    # Create solver
    solver = BDF2Solver(state, transport, bc_handler, params.get('time_step', {}))
    
    # Set Dirichlet values for H2, HeI
    if 'nH2_bc' in params:
        solver.set_dirichlet_value('H2', params['nH2_bc'])
    if 'nHeI_bc' in params:
        solver.set_dirichlet_value('HeI', params['nHeI_bc'])
    
    # Output callback
    output_times = []
    def output_callback(solver, t):
        output_interval = params.get('output_interval', 1e-4)
        if len(output_times) == 0 or t - output_times[-1] >= output_interval:
            output_times.append(t)
            if output_dir:
                write_csv_output(solver.state, t, radial_positions, output_dir)
            print(f"t = {t:.6e} s, dt = {solver.dt:.6e} s")
    
    # Run simulation
    t_end = params.get('tmainend', 1e-3)
    solver.run_until(t_end, callback=output_callback)
    
    return state
