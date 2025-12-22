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
        name: str = ""
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
        """
        self.V_space = V_space
        self.D_model = D_model
        self.V_model = V_model
        
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
        
        This is a placeholder for implementing Bohm/classical models.
        Currently only supports fixed coefficients.
        
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
            # TODO: Implement Bohm diffusion
            # D_Bohm = T / (16 * e * B) in SI units
            # D_Bohm [m²/s] = T[eV] / (16 * B[T])
            self._update_bohm_diffusion(species, B_field)
        elif self.D_model == DiffusionModel.CLASSICAL:
            # TODO: Implement classical transport
            self._update_classical_diffusion(state, species, B_field)
        
        if self.V_model == ConvectionModel.FIXED:
            # Keep current values
            pass
        elif self.V_model == ConvectionModel.PRESSURE:
            # TODO: Implement pressure-driven convection
            self._update_pressure_convection(state, species, B_field)
    
    def _update_bohm_diffusion(
        self,
        species: Optional[Species],
        B_field: Optional[np.ndarray]
    ) -> None:
        """
        Update diffusion using Bohm model: D = T / (16 * B).
        
        Placeholder implementation - to be completed.
        """
        if species is None or B_field is None:
            return
        
        T = species.T  # Temperature in eV
        # D_Bohm [m²/s] = T[eV] / (16 * B[T])
        D_bohm = T / (16.0 * B_field)
        
        # Apply with minimum value for stability
        D_min = 0.01  # Minimum diffusion [m²/s]
        self.D.x.array[:] = np.maximum(D_bohm, D_min)
    
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
        V_model: ConvectionModel = ConvectionModel.FIXED
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
            
        Returns
        -------
        coeff : TransportCoefficients
            Transport coefficient object for this species.
        """
        coeff = TransportCoefficients(
            self.V_space, D_value, V_value, D_model, V_model, name
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
        
        # Determine models from params
        diff_type = params.get("diffusion_type", "fixed").lower()
        conv_type = params.get("convection_type", "fixed").lower()
        
        D_model = DiffusionModel.FIXED
        if diff_type == "bohm":
            D_model = DiffusionModel.BOHM
        elif diff_type == "classical":
            D_model = DiffusionModel.CLASSICAL
        
        V_model = ConvectionModel.FIXED
        if conv_type == "pressure":
            V_model = ConvectionModel.PRESSURE
        
        # Initialize coefficients for each species
        for name in species_names:
            # Check for species-specific overrides
            D_species = params.get(f"D_{name}", D_fix)
            V_species = params.get(f"V_{name}", V_fix)
            
            self.add_species(name, D_species, V_species, D_model, V_model)
    
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


def compute_neutral_diffusion(
    n_background: np.ndarray,
    T_neutral: np.ndarray,
    m_neutral: float,
    sigma: float = 1e-19
) -> np.ndarray:
    """
    Compute neutral diffusion coefficient from collision frequency.
    
    D = v_th² / nu_coll = v_th * lambda_mfp
    
    Parameters
    ----------
    n_background : np.ndarray
        Background density for collisions [m^-3].
    T_neutral : np.ndarray
        Neutral temperature [eV].
    m_neutral : float
        Neutral mass [amu].
    sigma : float
        Collision cross-section [m²].
        
    Returns
    -------
    D : np.ndarray
        Diffusion coefficient [m²/s].
    """
    EV_TO_J = 1.60218e-19
    AMU_TO_KG = 1.66054e-27
    
    m_kg = m_neutral * AMU_TO_KG
    T_J = T_neutral * EV_TO_J
    
    # Thermal velocity
    v_th = np.sqrt(2 * T_J / m_kg)
    
    # Collision frequency
    nu = n_background * sigma * v_th
    
    # Diffusion coefficient (avoid division by zero)
    with np.errstate(divide='ignore', invalid='ignore'):
        D = np.where(nu > 1e-10, v_th**2 / nu, 1e10)
    
    # Limit maximum diffusion
    D = np.minimum(D, 1e4)
    
    return D
