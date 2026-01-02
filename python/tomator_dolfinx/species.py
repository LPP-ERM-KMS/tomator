"""
Species definitions for plasma transport simulation.

Defines the Species dataclass and PlasmaState container for managing
all plasma species (ions, neutrals, electrons) in the simulation.
"""

from dataclasses import dataclass, field
from typing import Optional, List, Dict
import numpy as np
from dolfinx import fem
from dolfinx.mesh import Mesh


# Physical constants
AMU_TO_KG = 1.66054e-27  # atomic mass unit to kg
EV_TO_J = 1.60218e-19    # eV to Joules
ME_KG = 9.10938e-31      # electron mass in kg


@dataclass
class Species:
    """
    Represents a plasma species with density and energy fields.
    
    Attributes
    ----------
    name : str
        Species identifier (e.g., "H", "Hi", "e").
    mass : float
        Mass in atomic mass units (amu). Electrons use 5.4858e-4.
    charge : int
        Charge state (0 for neutrals, 1 for singly ionized, etc., -1 for electrons).
    V : fem.FunctionSpace
        Function space for this species.
    n : fem.Function
        Particle density field [m^-3].
    E : fem.Function
        Energy density field [eV * m^-3] (E = 3/2 * n * T).
    n_prev : fem.Function
        Density at previous time step (for BDF2).
    n_prev2 : fem.Function
        Density at two time steps back (for BDF2).
    E_prev : fem.Function
        Energy at previous time step.
    E_prev2 : fem.Function
        Energy at two time steps back.
    solve_density : bool
        Whether to solve density equation (False for electrons, computed from quasi-neutrality).
    solve_energy : bool
        Whether to solve energy equation.
    bc_type : str
        Boundary condition type: "robin" (decay length) or "dirichlet".
    """
    name: str
    mass: float  # in amu
    charge: int
    V: fem.FunctionSpace
    n: fem.Function = field(init=False)
    E: fem.Function = field(init=False)
    n_prev: fem.Function = field(init=False)
    n_prev2: fem.Function = field(init=False)
    E_prev: fem.Function = field(init=False)
    E_prev2: fem.Function = field(init=False)
    solve_density: bool = True
    solve_energy: bool = True
    bc_type: str = "robin"  # "robin" or "dirichlet"
    
    def __post_init__(self):
        """Initialize function fields after dataclass creation."""
        self.n = fem.Function(self.V, name=f"n_{self.name}")
        self.E = fem.Function(self.V, name=f"E_{self.name}")
        self.n_prev = fem.Function(self.V, name=f"n_{self.name}_prev")
        self.n_prev2 = fem.Function(self.V, name=f"n_{self.name}_prev2")
        self.E_prev = fem.Function(self.V, name=f"E_{self.name}_prev")
        self.E_prev2 = fem.Function(self.V, name=f"E_{self.name}_prev2")
    
    @property
    def T(self) -> np.ndarray:
        """
        Temperature in eV, computed from E = 3/2 * n * T.
        
        Returns
        -------
        T : np.ndarray
            Temperature array [eV]. Clamped to [T_MIN, T_MAX] for stability.
        """
        # Temperature clamp limits (matching collisions.py)
        T_MIN = 0.026  # Room temperature ~0.026 eV
        T_MAX_HEAVY = 1000.0  # Max for ions/neutrals
        T_MAX_ELECTRON = 2e4  # Max for electrons
        
        n_arr = self.n.x.array
        E_arr = self.E.x.array
        # Avoid division by zero and ensure non-negative temperature
        with np.errstate(divide='ignore', invalid='ignore'):
            T = np.where(n_arr > 1e-10, E_arr / (1.5 * n_arr), T_MIN)
            # Clamp to physical range
            T_max = T_MAX_ELECTRON if self.name == "e" else T_MAX_HEAVY
            T = np.clip(T, T_MIN, T_max)
        return T
    
    @property
    def mass_kg(self) -> float:
        """Mass in kg."""
        return self.mass * AMU_TO_KG
    
    @property
    def thermal_velocity(self) -> np.ndarray:
        """
        Thermal velocity v_th = sqrt(2 * T / m) in m/s.
        
        Returns
        -------
        v_th : np.ndarray
            Thermal velocity array [m/s].
        """
        T_J = self.T * EV_TO_J  # Convert eV to Joules
        return np.sqrt(2 * T_J / self.mass_kg)
    
    def initialize_uniform(self, n0: float, T0: float) -> None:
        """
        Initialize with uniform density and temperature.
        
        Parameters
        ----------
        n0 : float
            Initial density [m^-3].
        T0 : float
            Initial temperature [eV].
        """
        E0 = 1.5 * n0 * T0
        self.n.x.array[:] = n0
        self.E.x.array[:] = E0
        self.n_prev.x.array[:] = n0
        self.n_prev2.x.array[:] = n0
        self.E_prev.x.array[:] = E0
        self.E_prev2.x.array[:] = E0
    
    def initialize_profile(self, n_profile: np.ndarray, T_profile: np.ndarray) -> None:
        """
        Initialize with given density and temperature profiles.
        
        Parameters
        ----------
        n_profile : np.ndarray
            Density profile [m^-3].
        T_profile : np.ndarray
            Temperature profile [eV].
        """
        E_profile = 1.5 * n_profile * T_profile
        self.n.x.array[:] = n_profile
        self.E.x.array[:] = E_profile
        self.n_prev.x.array[:] = n_profile
        self.n_prev2.x.array[:] = n_profile
        self.E_prev.x.array[:] = E_profile
        self.E_prev2.x.array[:] = E_profile
    
    def store_previous(self) -> None:
        """Store current solution for BDF2 time stepping."""
        self.n_prev2.x.array[:] = self.n_prev.x.array
        self.n_prev.x.array[:] = self.n.x.array
        self.E_prev2.x.array[:] = self.E_prev.x.array
        self.E_prev.x.array[:] = self.E.x.array
    
    def max_relative_change(self) -> float:
        """
        Compute maximum relative change from previous time step.
        
        Returns
        -------
        max_change : float
            Maximum of |n - n_prev| / n_prev over all mesh points.
        """
        n_arr = self.n.x.array
        n_prev_arr = self.n_prev.x.array
        
        with np.errstate(divide='ignore', invalid='ignore'):
            rel_change = np.abs(n_arr - n_prev_arr) / np.maximum(n_prev_arr, 1e-20)
        
        return np.max(rel_change)


class PlasmaState:
    """
    Container for all plasma species in the simulation.
    
    Manages electrons, hydrogen species (H, H+), and optionally
    helium (HeI, HeII, HeIII), molecular hydrogen (H2, H2+, H3+),
    and carbon impurities.
    
    Attributes
    ----------
    mesh : dolfinx.mesh.Mesh
        The computational mesh.
    V : fem.FunctionSpace
        Shared function space for all species.
    electrons : Species
        Electron species.
    species : Dict[str, Species]
        Dictionary of all species by name.
    """
    
    def __init__(self, mesh: Mesh, degree: int = 1):
        """
        Initialize plasma state with mesh and function space.
        
        Parameters
        ----------
        mesh : dolfinx.mesh.Mesh
            1D interval mesh.
        degree : int
            Polynomial degree for function space (default: 1 = linear).
        """
        self.mesh = mesh
        self.V = fem.functionspace(mesh, ("Lagrange", degree))
        self.species: Dict[str, Species] = {}
        
        # Initialize core species (electrons, H, H+)
        self._init_core_species()
    
    def _init_core_species(self) -> None:
        """Initialize electrons, atomic hydrogen, and H+ ions."""
        # Electrons (density from quasi-neutrality, solve energy only)
        self.electrons = Species(
            name="e",
            mass=5.4858e-4,  # electron mass in amu
            charge=-1,
            V=self.V,
            solve_density=False,  # Computed from quasi-neutrality
            solve_energy=True,
            bc_type="robin"
        )
        self.species["e"] = self.electrons
        
        # Atomic hydrogen (neutral)
        self.species["H"] = Species(
            name="H",
            mass=1.00784,
            charge=0,
            V=self.V,
            bc_type="robin"
        )
        
        # H+ (proton)
        self.species["Hi"] = Species(
            name="Hi",
            mass=1.00784,
            charge=1,
            V=self.V,
            bc_type="robin"
        )
    
    def add_helium_species(self) -> None:
        """Add helium species: HeI (neutral), HeII (He+), HeIII (He2+)."""
        # HeI - neutral helium (Dirichlet BC - fixed edge density)
        self.species["HeI"] = Species(
            name="HeI",
            mass=4.0026,
            charge=0,
            V=self.V,
            bc_type="dirichlet"
        )
        
        # HeII - singly ionized helium (Robin BC)
        self.species["HeII"] = Species(
            name="HeII",
            mass=4.0026,
            charge=1,
            V=self.V,
            bc_type="robin"
        )
        
        # HeIII - fully ionized helium (Robin BC)
        self.species["HeIII"] = Species(
            name="HeIII",
            mass=4.0026,
            charge=2,
            V=self.V,
            bc_type="robin"
        )
    
    def add_molecular_hydrogen_species(self) -> None:
        """Add molecular hydrogen species: H2, H2+, H3+."""
        # H2 - molecular hydrogen (Dirichlet BC)
        self.species["H2"] = Species(
            name="H2",
            mass=2.01568,
            charge=0,
            V=self.V,
            bc_type="dirichlet"
        )
        
        # H2+ ion (Robin BC)
        self.species["H2i"] = Species(
            name="H2i",
            mass=2.01568,
            charge=1,
            V=self.V,
            bc_type="robin"
        )
        
        # H3+ ion (Robin BC)
        self.species["H3i"] = Species(
            name="H3i",
            mass=3.02352,
            charge=1,
            V=self.V,
            bc_type="robin"
        )
    
    def add_carbon_species(self) -> None:
        """Add carbon impurity species: CI, CII, CIII, CIV."""
        carbon_mass = 12.011
        
        # CI - neutral carbon
        self.species["CI"] = Species(
            name="CI",
            mass=carbon_mass,
            charge=0,
            V=self.V,
            bc_type="robin"
        )
        
        # CII through CIV (ionized states)
        for i, charge in enumerate([1, 2, 3], start=2):
            name = f"C{['I', 'II', 'III', 'IV'][i-1]}"
            self.species[name] = Species(
                name=name,
                mass=carbon_mass,
                charge=charge,
                V=self.V,
                bc_type="robin"
            )
    
    @property
    def ions(self) -> List[Species]:
        """List of all ion species (charge > 0)."""
        return [s for s in self.species.values() if s.charge > 0]
    
    @property
    def neutrals(self) -> List[Species]:
        """List of all neutral species (charge = 0)."""
        return [s for s in self.species.values() if s.charge == 0]
    
    @property
    def all_species(self) -> List[Species]:
        """List of all species."""
        return list(self.species.values())
    
    def compute_electron_density(self) -> None:
        """
        Compute electron density from quasi-neutrality.
        
        n_e = sum(Z_i * n_i) for all ion species.
        """
        ne = np.zeros_like(self.electrons.n.x.array)
        
        for species in self.ions:
            ne += species.charge * species.n.x.array
        
        # Ensure non-negative
        ne = np.maximum(ne, 1e-10)
        self.electrons.n.x.array[:] = ne
    
    def initialize_from_params(self, params: dict) -> None:
        """
        Initialize all species from parameter dictionary.
        
        Uses Gaussian spatial profiles matching C++ Tomator1D:
        - fct_n0 = background + exp(-((r - r_max) / width)^2)
        
        Parameters
        ----------
        params : dict
            Dictionary with initial conditions including profile parameters.
        """
        # Get profile parameters (converted to SI in json_input.py)
        rmaxini = params.get("rmaxini", 0.92)  # Peak location [m]
        widthini = params.get("widthini", 0.05)  # Gaussian width [m]
        nebackgroundl = params.get("nebackgroundl", 1e-5)  # HFS background (fraction)
        nebackgroundr = params.get("nebackgroundr", 1e-3)  # LFS background (fraction)
        
        # Get mesh coordinates
        coords = self.mesh.geometry.x[:, 0]  # 1D coordinates [m]
        
        # Compute Gaussian profile (matching C++ functions.cpp)
        # fct_n0 = background + exp(-((r - rmaxini)/widthini)^2)
        # Different widths for HFS (r < rmaxini) and LFS (r > rmaxini)
        fct_n0 = np.zeros_like(coords)
        for i, r in enumerate(coords):
            if r <= rmaxini:
                fct_n0[i] = nebackgroundl + np.exp(-((r - rmaxini) / widthini)**2)
            else:
                # C++ uses 2*width on LFS side
                fct_n0[i] = nebackgroundr + np.exp(-((r - rmaxini) / (2 * widthini))**2)
        
        # Temperatures [eV]
        Te0 = params.get("Te0", 10.0)
        Ta0 = params.get("Ta0", 0.026)  # Neutral temp / room temp
        
        # Initialize ions with Gaussian profiles
        # H+ (proton)
        if "Hi" in self.species:
            nHi0 = params.get("nHi0", 1e18)  # Peak density [m^-3]
            THi0 = params.get("THi0", Te0)
            n_profile = nHi0 * fct_n0
            T_profile = np.full_like(coords, THi0)
            self.species["Hi"].initialize_profile(n_profile, T_profile)
        
        # HeII - He+ (Gaussian profile)
        if "HeII" in self.species:
            nHeII0 = params.get("nHeII0", 1e18)
            THeII0 = params.get("THeII0", Te0)
            n_profile = nHeII0 * fct_n0
            T_profile = np.full_like(coords, THeII0)
            self.species["HeII"].initialize_profile(n_profile, T_profile)
        
        # HeIII - He++ (Gaussian profile)
        if "HeIII" in self.species:
            nHeIII0 = params.get("nHeIII0", 1e13)
            THeIII0 = params.get("THeIII0", Te0)
            n_profile = nHeIII0 * fct_n0
            T_profile = np.full_like(coords, THeIII0)
            self.species["HeIII"].initialize_profile(n_profile, T_profile)
        
        # Molecular ions with Gaussian profiles
        if "H2i" in self.species:
            nH2i0 = params.get("nH2i0", 1e12)
            TH2i0 = params.get("TH2i0", Te0)
            n_profile = nH2i0 * fct_n0
            T_profile = np.full_like(coords, TH2i0)
            self.species["H2i"].initialize_profile(n_profile, T_profile)
        
        if "H3i" in self.species:
            nH3i0 = params.get("nH3i0", 1e12)
            TH3i0 = params.get("TH3i0", Te0)
            n_profile = nH3i0 * fct_n0
            T_profile = np.full_like(coords, TH3i0)
            self.species["H3i"].initialize_profile(n_profile, T_profile)
        
        # Neutrals: UNIFORM (not Gaussian) - matches C++ which uses nH0 * fct_n0 for H
        # but sets nH2[im] = nH20 and nHeI[im] = nHeI0 (uniform)
        # Actually, C++ uses fct_n0 for H but uniform for H2 and HeI
        if "H" in self.species:
            nH0 = params.get("nH0", 1e16)
            TH0 = params.get("TH0", Ta0)
            n_profile = nH0 * fct_n0  # H also uses Gaussian in C++
            T_profile = np.full_like(coords, TH0)
            self.species["H"].initialize_profile(n_profile, T_profile)
        
        if "H2" in self.species:
            nH20 = params.get("nH20", 1e15)
            TH20 = params.get("TH20", Ta0)
            self.species["H2"].initialize_uniform(nH20, TH20)  # Uniform
        
        if "HeI" in self.species:
            nHeI0 = params.get("nHeI0", 1e15)
            THeI0 = params.get("THeI0", Ta0)
            self.species["HeI"].initialize_uniform(nHeI0, THeI0)  # Uniform
        
        # Compute electron density from quasi-neutrality
        self.compute_electron_density()
        
        # Set electron energy from given temperature (uniform Te)
        ne_arr = self.electrons.n.x.array
        self.electrons.E.x.array[:] = 1.5 * ne_arr * Te0
        self.electrons.E_prev.x.array[:] = self.electrons.E.x.array
        self.electrons.E_prev2.x.array[:] = self.electrons.E.x.array
    
    def store_all_previous(self) -> None:
        """Store current solution for all species (for BDF2)."""
        for species in self.all_species:
            species.store_previous()
    
    def max_relative_change(self) -> float:
        """
        Get maximum relative change across all species.
        
        Returns
        -------
        max_change : float
            Maximum relative change for timestep adaptation.
        """
        changes = [s.max_relative_change() for s in self.all_species if s.solve_density]
        return max(changes) if changes else 0.0
