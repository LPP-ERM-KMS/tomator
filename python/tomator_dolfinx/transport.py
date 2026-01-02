"""
Transport coefficient calculations for plasma species.

Provides diffusion (D) and convection (V) coefficients as dolfinx Functions,
supporting uniform values initially with placeholders for Bohm/classical models.
"""

from typing import Optional, Dict, Tuple
from enum import Enum
import numpy as np
from dolfinx import fem, default_scalar_type
from dolfinx.mesh import Mesh

from .species import Species, PlasmaState


class DiffusionModel(Enum):
    """Available diffusion coefficient models."""
    FIXED = "fixed"         # Spatially uniform, fixed value
    BOHM = "bohm"           # D ~ T / B
    CLASSICAL = "classical" # Classical transport with collisions


class ConvectionModel(Enum):
    """Available convection velocity models."""
    FIXED = "fixed"         # Spatially uniform, fixed value
    PRESSURE = "pressure"   # Pressure-driven pinch


class TransportCoefficients:
    """
    Transport coefficients (D, V) for a species.
    
    Stores diffusion and convection as dolfinx Functions that can be
    spatially varying. Initially supports uniform values; models for
    Bohm, classical transport to be added later.
    
    Attributes
    ----------
    D : fem.Function
        Diffusion coefficient [m²/s].
    V : fem.Function
        Convection velocity [m/s] (positive = outward).
    D_model : DiffusionModel
        Model used for diffusion coefficient.
    V_model : ConvectionModel
        Model used for convection velocity.
    """
    
    def __init__(
        self,
        V_space: fem.FunctionSpace,
        D_value: float = 1.0,
        V_value: float = 0.0,
        D_model: DiffusionModel = DiffusionModel.FIXED,
        V_model: ConvectionModel = ConvectionModel.FIXED,
        name: str = "",
        Dfsave: float = 1.0
    ):
        """
        Initialize transport coefficients.
        
        Parameters
        ----------
        V_space : fem.FunctionSpace
            Function space for coefficient fields.
        D_value : float
            Initial/fixed diffusion coefficient [m²/s].
        V_value : float
            Initial/fixed convection velocity [m/s].
        D_model : DiffusionModel
            Diffusion model type.
        V_model : ConvectionModel
            Convection model type.
        name : str
            Species name for labeling functions.
        Dfsave : float
            Bohm diffusion scaling factor (default 1.0).
        """
        self.V_space = V_space
        self.D_model = D_model
        self.V_model = V_model
        self.Dfsave = Dfsave
        
        # Create coefficient Functions
        self.D = fem.Function(V_space, name=f"D_{name}" if name else "D")
        self.V = fem.Function(V_space, name=f"V_{name}" if name else "V")
        
        # Store base values for models
        self._D_base = D_value
        self._V_base = V_value
        
        # Initialize with uniform values
        self.D.x.array[:] = D_value
        self.V.x.array[:] = V_value
    
    def set_uniform(self, D_value: float, V_value: float) -> None:
        """
        Set uniform (spatially constant) coefficients.
        
        Parameters
        ----------
        D_value : float
            Diffusion coefficient [m²/s].
        V_value : float
            Convection velocity [m/s].
        """
        self.D.x.array[:] = D_value
        self.V.x.array[:] = V_value
        self._D_base = D_value
        self._V_base = V_value
    
    def update(
        self,
        state: Optional[PlasmaState] = None,
        species: Optional[Species] = None,
        B_field: Optional[np.ndarray] = None
    ) -> None:
        """
        Update coefficients based on current plasma state.
        
        Parameters
        ----------
        state : PlasmaState, optional
            Current plasma state (for temperature-dependent models).
        species : Species, optional
            Species for which to compute coefficients.
        B_field : np.ndarray, optional
            Magnetic field strength array [T].
        """
        if self.D_model == DiffusionModel.FIXED:
            # Keep current values
            pass
        elif self.D_model == DiffusionModel.BOHM:
            # Bohm diffusion: D = Dfsave * T_eff / B
            self._update_bohm_diffusion(state, B_field, self.Dfsave)
        elif self.D_model == DiffusionModel.CLASSICAL:
            # Classical transport with collisions
            self._update_classical_diffusion(state, species, B_field)
        
        if self.V_model == ConvectionModel.FIXED:
            # Keep current values
            pass
        elif self.V_model == ConvectionModel.PRESSURE:
            # Pressure-driven convection
            self._update_pressure_convection(state, species, B_field)
    
    def _update_bohm_diffusion(
        self,
        state,  # PlasmaState
        B_field: Optional[np.ndarray],
        Dfsave: float = 1.0
    ) -> None:
        """
        Update diffusion using Bohm model: D = Dfsave * T_eff / B.
        
        Uses pressure-weighted average temperature from all ion species,
        matching the C++ Tomator1D implementation.
        
        Parameters
        ----------
        state : PlasmaState
            Current plasma state with densities and temperatures.
        B_field : np.ndarray
            Radial magnetic field [T].
        Dfsave : float
            Diffusion scaling factor.
        """
        if state is None or B_field is None:
            return
        
        D_bohm = compute_bohm_diffusion_from_state(
            state, B_field, Dfsave=Dfsave,
            D_min=0.01, D_max=1e4
        )
        self.D.x.array[:] = D_bohm
    
    def _update_classical_diffusion(
        self,
        state: Optional[PlasmaState],
        species: Optional[Species],
        B_field: Optional[np.ndarray]
    ) -> None:
        """
        Update diffusion using classical transport model.
        
        Placeholder implementation - to be completed.
        Classical: D = (1/3) * nu * lambda_mfp * (rho_s + lambda_mfp * B_h/B_t)
        """
        # TODO: Implement classical transport from C++ transport.cpp
        pass
    
    def _update_pressure_convection(
        self,
        state: Optional[PlasmaState],
        species: Optional[Species],
        B_field: Optional[np.ndarray]
    ) -> None:
        """
        Update convection using pressure-driven model.
        
        Placeholder implementation - to be completed.
        """
        # TODO: Implement pressure-driven pinch from C++ transport.cpp
        pass


class TransportManager:
    """
    Manages transport coefficients for all species in the simulation.
    
    Provides unified interface for initializing and updating transport
    coefficients for ions, neutrals, and electrons.
    
    Attributes
    ----------
    coefficients : Dict[str, TransportCoefficients]
        Transport coefficients keyed by species name.
    """
    
    def __init__(self, V_space: fem.FunctionSpace):
        """
        Initialize transport manager.
        
        Parameters
        ----------
        V_space : fem.FunctionSpace
            Function space for coefficient fields.
        """
        self.V_space = V_space
        self.coefficients: Dict[str, TransportCoefficients] = {}
    
    def add_species(
        self,
        name: str,
        D_value: float = 1.0,
        V_value: float = 0.0,
        D_model: DiffusionModel = DiffusionModel.FIXED,
        V_model: ConvectionModel = ConvectionModel.FIXED,
        Dfsave: float = 1.0
    ) -> TransportCoefficients:
        """
        Add transport coefficients for a species.
        
        Parameters
        ----------
        name : str
            Species name.
        D_value : float
            Initial diffusion coefficient [m²/s].
        V_value : float
            Initial convection velocity [m/s].
        D_model : DiffusionModel
            Diffusion model type.
        V_model : ConvectionModel
            Convection model type.
        Dfsave : float
            Bohm diffusion scaling factor (default 1.0).
            
        Returns
        -------
        coeff : TransportCoefficients
            Transport coefficient object for this species.
        """
        coeff = TransportCoefficients(
            self.V_space, D_value, V_value, D_model, V_model, name, Dfsave
        )
        self.coefficients[name] = coeff
        return coeff
    
    def get(self, name: str) -> TransportCoefficients:
        """
        Get transport coefficients for a species.
        
        Parameters
        ----------
        name : str
            Species name.
            
        Returns
        -------
        coeff : TransportCoefficients
            Transport coefficients for the species.
        """
        return self.coefficients[name]
    
    def initialize_from_params(
        self,
        params: dict,
        species_names: list
    ) -> None:
        """
        Initialize transport coefficients from parameter dictionary.
        
        Parameters
        ----------
        params : dict
            Dictionary with transport parameters, e.g.:
            {
                "Dfix": 1.0,
                "Vfix": 0.0,
                "diffusion_type": "fixed",
                "convection_type": "fixed",
                ...
            }
        species_names : list
            List of species names to initialize.
        """
        # Get global values (apply to all species unless overridden)
        D_fix = params.get("Dfix", 1.0)
        V_fix = params.get("Vfix", 0.0)
        Dfsave = params.get("Dfsave", params.get("Dfact", 1.0))  # Bohm scaling factor
        
        # Determine models from params - check C++ style bool flags first
        # Priority: bDfix > bDbohm > bDscaling > diffusion_type string
        D_model = DiffusionModel.FIXED  # Default
        
        if params.get("bDfix", True):
            # Fixed diffusion coefficient (default behavior)
            D_model = DiffusionModel.FIXED
        elif params.get("bDbohm", False):
            # Bohm diffusion: D = Dfsave * T_eff / B
            D_model = DiffusionModel.BOHM
        elif params.get("bDscaling", False):
            # Classical diffusion with collision scaling
            D_model = DiffusionModel.CLASSICAL
        else:
            # Fallback to string-based type if no bool flags
            diff_type = params.get("diffusion_type", "fixed").lower()
            if diff_type == "bohm":
                D_model = DiffusionModel.BOHM
            elif diff_type == "classical":
                D_model = DiffusionModel.CLASSICAL
        
        # Convection model
        V_model = ConvectionModel.FIXED
        conv_type = params.get("convection_type", "fixed").lower()
        if conv_type == "pressure":
            V_model = ConvectionModel.PRESSURE
        
        # Initialize coefficients for each species
        for name in species_names:
            # Check for species-specific overrides
            D_species = params.get(f"D_{name}", D_fix)
            V_species = params.get(f"V_{name}", V_fix)
            
            self.add_species(name, D_species, V_species, D_model, V_model, Dfsave)
    
    def update_all(
        self,
        state: PlasmaState,
        B_field: Optional[np.ndarray] = None
    ) -> None:
        """
        Update all transport coefficients based on current state.
        
        Parameters
        ----------
        state : PlasmaState
            Current plasma state.
        B_field : np.ndarray, optional
            Magnetic field strength array [T].
        """
        for name, coeff in self.coefficients.items():
            species = state.species.get(name)
            coeff.update(state, species, B_field)


# =========================================================================
# Physics-Based Diffusion from C++ Tomator1D (transport.cpp)
# 
# D = (1/3) * vth / (1/mfp + 1/(a/2))
#
# This uses a harmonic mean of collision-limited (mfp) and geometry-limited (a/2)
# mean free paths. At high collisionality, D ~ vth * mfp / 3.
# At low collisionality, D is limited by geometry to ~ vth * a/6.
#
# vth = sqrt(8 * kB * T / (pi * m))  [mean thermal speed]
# mfp = vth / nu                     [collision mean free path]
# nu = sum of all collision frequencies [1/s]
# =========================================================================

# Physical constants
_EV_TO_J = 1.60218e-19
_AMU_TO_KG = 1.66054e-27


def compute_thermal_velocity(T_eV: np.ndarray, m_amu: float) -> np.ndarray:
    """
    Compute mean thermal velocity.
    
    vth = sqrt(8 * kB * T / (pi * m))
    
    Parameters
    ----------
    T_eV : np.ndarray
        Temperature [eV].
    m_amu : float
        Mass [amu].
        
    Returns
    -------
    vth : np.ndarray
        Mean thermal velocity [m/s].
    """
    m_kg = m_amu * _AMU_TO_KG
    T_J = np.maximum(T_eV, 0.01) * _EV_TO_J  # Minimum 0.01 eV
    
    return np.sqrt(8.0 * T_J / (np.pi * m_kg))


def compute_neutral_diffusion(
    nu_collision: np.ndarray,
    T_eV: np.ndarray,
    m_amu: float,
    a_minor: float,
    D_min: float = 1e-4,
    D_max: float = 1e4
) -> np.ndarray:
    """
    Compute neutral diffusion coefficient from collision frequency.
    
    From C++ Tomator1D (transport.cpp):
    D = (1/3) * vth / (1/mfp + 1/(a/2))
    
    This is a harmonic mean of:
    - Collisional limit: D ~ vth * mfp / 3 = vth² / (3 * nu)
    - Geometric limit: D ~ vth * a / 6
    
    Parameters
    ----------
    nu_collision : np.ndarray
        Total collision frequency [1/s].
    T_eV : np.ndarray
        Species temperature [eV].
    m_amu : float
        Species mass [amu].
    a_minor : float
        Minor radius (characteristic length scale) [m].
    D_min : float
        Minimum diffusion coefficient [m²/s].
    D_max : float
        Maximum diffusion coefficient [m²/s].
        
    Returns
    -------
    D : np.ndarray
        Diffusion coefficient [m²/s].
    """
    # Thermal velocity
    vth = compute_thermal_velocity(T_eV, m_amu)
    
    # Mean free path: mfp = vth / nu
    # Avoid division by zero
    nu_safe = np.maximum(nu_collision, 1e-10)
    mfp = vth / nu_safe
    
    # Geometry limit: a/2
    L_geom = 0.5 * a_minor
    
    # Harmonic mean: 1/L_eff = 1/mfp + 1/L_geom
    # D = (1/3) * vth / (1/mfp + 1/L_geom) = (1/3) * vth * L_eff
    # where L_eff = 1 / (1/mfp + 1/L_geom) = mfp * L_geom / (mfp + L_geom)
    L_eff = mfp * L_geom / (mfp + L_geom)
    
    D = (1.0/3.0) * vth * L_eff
    
    # Apply bounds
    D = np.clip(D, D_min, D_max)
    
    return D


def compute_neutral_diffusion_from_state(
    species_name: str,
    state,  # PlasmaState
    a_minor: float,
    include_self_collision: bool = True,
    nu_collision: dict = None
) -> np.ndarray:
    """
    Compute physics-based diffusion for neutral species from plasma state.
    
    This function is for NEUTRALS ONLY (H, H2, HeI). Ion diffusion uses
    different physics (Bohm or fixed coefficients).
    
    Convenience function that:
    1. Extracts densities and temperatures from PlasmaState
    2. Uses collision frequency for the species (from nu_collision dict)
    3. Computes diffusion coefficient using harmonic mean of mfp and geometry
    
    Parameters
    ----------
    species_name : str
        Neutral species name ('H', 'H2', 'HeI').
    state : PlasmaState
        Current plasma state with all densities and temperatures.
    a_minor : float
        Minor radius [m].
    include_self_collision : bool
        Whether to include self-collisions.
    nu_collision : dict, optional
        Pre-computed collision frequencies from compute_collision_sources.
        If provided, uses nu[species_name] instead of computing locally.
        
    Returns
    -------
    D : np.ndarray
        Diffusion coefficient [m²/s].
        
    Raises
    ------
    ValueError
        If species_name is not a neutral species (H, H2, HeI).
    """
    # Only neutrals use this physics-based diffusion
    neutral_species = {'H', 'H2', 'HeI'}
    if species_name not in neutral_species:
        raise ValueError(f"compute_neutral_diffusion_from_state is for neutrals only (H, H2, HeI), "
                         f"not '{species_name}'. Use Bohm or fixed diffusion for ions.")
    
    # Get species data
    def get_species_arrays(name):
        """Get density [m^-3] and temperature [eV] arrays."""
        s = state.species.get(name)
        if s is None:
            return None, None
        n = s.n.x.array
        T = s.T  # This is now always an ndarray from the Species.T property
        # Ensure T is an array of the same size as n
        if isinstance(T, np.ndarray):
            return n, T
        else:
            # Scalar fallback
            return n, np.full_like(n, float(T))
    
    # Temperature map for neutral species only
    _, TH = get_species_arrays('H')
    _, TH2 = get_species_arrays('H2')
    _, THeI = get_species_arrays('HeI')
    
    # Map species name to temperature
    T_map = {'H': TH, 'H2': TH2, 'HeI': THeI}
    T_species = T_map.get(species_name)
    if T_species is None:
        raise ValueError(f"Unknown neutral species: {species_name}")
    
    # Mass mapping for neutrals
    mass_map = {'H': 1.0, 'H2': 2.0, 'HeI': 4.0}
    m_amu = mass_map.get(species_name, 1.0)
    
    # Get collision frequency from pre-computed dict
    # nu_collision keys: 'e', 'H', 'Hi', 'H2', 'H2i', 'H3i', 'HeI', 'HeII', 'HeIII'
    if nu_collision is None or species_name not in nu_collision:
        # Fallback: use a small default nu to avoid zero division
        n_points = len(state.electrons.n.x.array)
        nu = np.full(n_points, 1e-10)
    else:
        nu = nu_collision[species_name]
    
    # Compute diffusion using physics-based model
    D = compute_neutral_diffusion(nu, T_species, m_amu, a_minor)
    
    return D


# =========================================================================
# Bohm Diffusion for Ions (from C++ transport.cpp)
#
# D_Bohm = Dfsave * T_eff / B_r
#
# where T_eff = (sum_i n_i * T_i) / n_e is the pressure-weighted average
# temperature of all ion species, and B_r is the radial magnetic field.
#
# Standard Bohm diffusion is D = T / (16 * e * B) in SI, or equivalently
# D [m²/s] = T[eV] / (16 * B[T]) ≈ 0.0625 * T/B
#
# The C++ code uses: Dionh = max(pressure_sum, ne*0.05) / ne
#                    Dion = max(1e2, Dfsave * Dionh / Br)
# where units are in [cm²/s] hence the 1e2 minimum (= 0.01 m²/s)
# =========================================================================


def compute_bohm_diffusion(
    ne: np.ndarray,
    Te: np.ndarray,
    n_ions: dict,
    T_ions: dict,
    B_radial: np.ndarray,
    Dfsave: float = 1.0,
    D_min: float = 0.01,
    D_max: float = 1e4,
    T_floor: float = 0.05
) -> np.ndarray:
    """
    Compute Bohm diffusion coefficient for ions.
    
    From C++ Tomator1D (transport.cpp), bDbohm mode:
    
    Dionh = max(ne*Te + sum(n_i * T_i), ne * T_floor) / ne
    D = max(D_min, Dfsave * Dionh / B_r)
    
    This is effectively D ~ T_eff / B where T_eff is the density-weighted
    average temperature.
    
    Parameters
    ----------
    ne : np.ndarray
        Electron density [m^-3].
    Te : np.ndarray
        Electron temperature [eV].
    n_ions : dict
        Dictionary of ion densities [m^-3] keyed by species name.
        Expected keys: 'Hi', 'H2i', 'H3i', 'HeII', 'HeIII'
    T_ions : dict
        Dictionary of ion temperatures [eV] keyed by species name.
        Expected keys: 'Hi', 'H2i', 'H3i', 'HeII', 'HeIII'
    B_radial : np.ndarray
        Radial magnetic field strength [T].
    Dfsave : float
        Diffusion scaling factor (default 1.0).
    D_min : float
        Minimum diffusion coefficient [m²/s] (default 0.01 = 1e2 cm²/s).
    D_max : float
        Maximum diffusion coefficient [m²/s].
    T_floor : float
        Minimum effective temperature [eV] (default 0.05).
        
    Returns
    -------
    D : np.ndarray
        Bohm diffusion coefficient [m²/s].
        
    Notes
    -----
    The C++ implementation uses units of [cm²/s] internally and has a minimum
    of 1e2 cm²/s = 0.01 m²/s. The Dfsave factor typically comes from tuning
    or PI control in the simulation.
    
    Standard Bohm diffusion is D = T / (16*B) [m²/s] with T in eV and B in T.
    The factor 1/16 comes from: D_Bohm = k_B T / (e B * 16)
    """
    # Compute pressure-weighted sum: sum(n_i * T_i) for all ion species
    pressure_sum = ne * Te  # Electron contribution
    
    # Add ion contributions
    ion_species = ['Hi', 'H2i', 'H3i', 'HeII', 'HeIII']
    for ion in ion_species:
        if ion in n_ions and ion in T_ions:
            n_i = n_ions[ion]
            T_i = T_ions[ion]
            # Ensure arrays are compatible
            if n_i is not None and T_i is not None:
                # Handle scalar or array temperatures
                if np.isscalar(T_i):
                    T_i = np.full_like(n_i, T_i)
                pressure_sum = pressure_sum + n_i * T_i
    
    # Effective temperature with floor: T_eff = max(pressure_sum, ne * T_floor) / ne
    ne_safe = np.maximum(ne, 1e10)  # Avoid division by zero
    T_eff = np.maximum(pressure_sum, ne_safe * T_floor) / ne_safe
    
    # Bohm diffusion: D = Dfsave * T_eff / B
    B_safe = np.maximum(np.abs(B_radial), 1e-6)  # Avoid division by zero
    D_bohm = Dfsave * T_eff / B_safe
    
    # Apply bounds
    D = np.clip(D_bohm, D_min, D_max)
    
    return D


def compute_bohm_diffusion_from_state(
    state,  # PlasmaState
    B_radial: np.ndarray,
    Dfsave: float = 1.0,
    D_min: float = 0.01,
    D_max: float = 1e4
) -> np.ndarray:
    """
    Compute Bohm diffusion coefficient from PlasmaState.
    
    Convenience function that extracts densities and temperatures from
    the PlasmaState object and calls compute_bohm_diffusion.
    
    Parameters
    ----------
    state : PlasmaState
        Current plasma state with all densities and temperatures.
    B_radial : np.ndarray
        Radial magnetic field strength [T].
    Dfsave : float
        Diffusion scaling factor (default 1.0).
    D_min : float
        Minimum diffusion coefficient [m²/s].
    D_max : float
        Maximum diffusion coefficient [m²/s].
        
    Returns
    -------
    D : np.ndarray
        Bohm diffusion coefficient [m²/s].
    """
    # Get electron density and temperature
    ne = state.electrons.n.x.array.copy()
    Te = state.electrons.T  # T property returns array
    
    # Build ion density and temperature dictionaries
    n_ions = {}
    T_ions = {}
    
    ion_species = ['Hi', 'H2i', 'H3i', 'HeII', 'HeIII']
    for ion in ion_species:
        s = state.species.get(ion)
        if s is not None:
            n_ions[ion] = s.n.x.array.copy()
            T_ions[ion] = s.T  # Property returns array
    
    return compute_bohm_diffusion(
        ne, Te, n_ions, T_ions, B_radial,
        Dfsave=Dfsave, D_min=D_min, D_max=D_max
    )
