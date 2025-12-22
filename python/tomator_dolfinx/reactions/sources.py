"""
Source term calculations for plasma species.

Computes collision source terms dn/dt and dE/dt for each species
based on the reaction rate coefficients.

NOTE: All densities are expected in SI units (m^-3).
Rate coefficients from rates.py are in CGS (cm³/s) and are converted internally.
Output source terms are in SI: dn [m^-3/s], dE [eV·m^-3/s].
"""

from typing import Dict, Tuple
import numpy as np

from .rates import ReactionRates


# Unit conversion: cm³/s to m³/s
CM3_TO_M3 = 1e-6
# Unit conversion: m^-3 to cm^-3
M3_TO_CM3 = 1e-6

# Maximum rate value to prevent overflow (1e30 is safe for float64)
MAX_RATE = 1e30


def ndamp(n_cgs: np.ndarray, nevac: float = 1.0) -> np.ndarray:
    """
    Density damping factor from C++ Tomator1D.
    
    ndamp(n) = (1 + nevac/n)^(-0.66)
    
    This limits reactions at very low densities (approaching vacuum).
    
    Parameters
    ----------
    n_cgs : np.ndarray
        Density [cm^-3].
    nevac : float
        Vacuum density reference [cm^-3]. Default is 1.0 from C++.
        
    Returns
    -------
    factor : np.ndarray
        Damping factor (0 to 1).
    """
    n_safe = np.maximum(n_cgs, 1e-10)  # Prevent division by zero
    return np.power(1.0 + nevac / n_safe, -0.66)


def safe_rate(k, *densities):
    """
    Compute reaction rate k * n1 * n2 * ... with overflow protection.
    
    Multiplies in order: k * n1, then result * n2, etc.
    Uses np.multiply with 'where' to handle potential overflow.
    Clips final result to prevent overflow.
    """
    with np.errstate(over='ignore', invalid='ignore'):
        result = np.asarray(k, dtype=np.float64)
        for n in densities:
            result = result * np.asarray(n, dtype=np.float64)
        # Replace inf/nan with max value
        result = np.where(np.isfinite(result), result, np.sign(result) * MAX_RATE)
        return np.clip(result, -MAX_RATE, MAX_RATE)


class SourceTerms:
    """
    Container for source term arrays for all species.
    
    Attributes
    ----------
    dn : Dict[str, np.ndarray]
        Density source terms [cm^-3/s] keyed by species name.
    dE : Dict[str, np.ndarray]
        Energy source terms [eV*cm^-3/s] keyed by species name.
    """
    
    def __init__(self, nmesh: int, species_names: list):
        """
        Initialize source term arrays.
        
        Parameters
        ----------
        nmesh : int
            Number of mesh points.
        species_names : list
            List of species names.
        """
        self.dn = {name: np.zeros(nmesh) for name in species_names}
        self.dE = {name: np.zeros(nmesh) for name in species_names}
        self.nmesh = nmesh
        self.species_names = species_names
    
    def reset(self) -> None:
        """Reset all source terms to zero."""
        for name in self.species_names:
            self.dn[name][:] = 0.0
            self.dE[name][:] = 0.0


def compute_collision_sources(
    ne: np.ndarray,
    Te: np.ndarray,
    nH: np.ndarray,
    TH: np.ndarray,
    nHi: np.ndarray,
    THi: np.ndarray,
    nHeI: np.ndarray = None,
    THeI: np.ndarray = None,
    nHeII: np.ndarray = None,
    THeII: np.ndarray = None,
    nHeIII: np.ndarray = None,
    THeIII: np.ndarray = None,
    use_ADAS: bool = True,
    include_H: bool = True,
    include_He: bool = True,
    include_cx: bool = True
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """
    Compute all collision source terms for density and energy.
    
    This is the main function that computes dn/dt and dE/dt for each species
    from all relevant reactions.
    
    Parameters
    ----------
    ne : np.ndarray
        Electron density [m^-3].
    Te : np.ndarray
        Electron temperature [eV].
    nH : np.ndarray
        Atomic hydrogen density [m^-3].
    TH : np.ndarray
        Atomic hydrogen temperature [eV].
    nHi : np.ndarray
        H+ density [m^-3].
    THi : np.ndarray
        H+ temperature [eV].
    nHeI : np.ndarray, optional
        Neutral helium density [m^-3].
    THeI : np.ndarray, optional
        Neutral helium temperature [eV].
    nHeII : np.ndarray, optional
        He+ density [m^-3].
    THeII : np.ndarray, optional
        He+ temperature [eV].
    nHeIII : np.ndarray, optional
        He++ density [m^-3].
    THeIII : np.ndarray, optional
        He++ temperature [eV].
    use_ADAS : bool
        Use ADAS cooling rates for He (True) or IAEA formulas (False).
    include_H : bool
        Include hydrogen reactions.
    include_He : bool
        Include helium reactions.
    include_cx : bool
        Include charge exchange reactions.
        
    Returns
    -------
    dn : Dict[str, np.ndarray]
        Density source terms [m^-3/s].
    dE : Dict[str, np.ndarray]
        Energy source terms [eV·m^-3/s].
    """
    nmesh = len(ne)
    
    # Initialize source arrays
    dn = {
        'e': np.zeros(nmesh),
        'H': np.zeros(nmesh),
        'Hi': np.zeros(nmesh),
    }
    dE = {
        'e': np.zeros(nmesh),
        'H': np.zeros(nmesh),
        'Hi': np.zeros(nmesh),
    }
    
    # Add He species if provided
    if nHeI is not None:
        dn['HeI'] = np.zeros(nmesh)
        dn['HeII'] = np.zeros(nmesh)
        dn['HeIII'] = np.zeros(nmesh)
        dE['HeI'] = np.zeros(nmesh)
        dE['HeII'] = np.zeros(nmesh)
        dE['HeIII'] = np.zeros(nmesh)
    
    rates = ReactionRates()
    
    # Convert SI densities to CGS for rate calculations (m^-3 -> cm^-3)
    ne_cgs = ne * M3_TO_CM3
    nH_cgs = nH * M3_TO_CM3
    nHi_cgs = nHi * M3_TO_CM3
    
    # -------------------------------------------------------------------
    # Electron-Hydrogen Reactions (bH)
    # Rate coefficients from rates.py are in cm³/s, convert to m³/s
    # ndamp() takes CGS densities
    # -------------------------------------------------------------------
    if include_H:
        # 1. H excitation: e + H -> e + H* (energy loss only)
        k_exc = rates.H_excitation(Te) * CM3_TO_M3
        rate_exc = safe_rate(k_exc, ne, nH)
        dE['e'] -= rate_exc * 10.2  # Excitation energy loss
        
        # 2. H ionization: e + H -> e + H+ + e
        k_ion = rates.H_ionization(Te) * ndamp(nH_cgs) * CM3_TO_M3
        rate_ion = safe_rate(k_ion, ne, nH)
        
        dn['H'] -= rate_ion
        dn['Hi'] += rate_ion
        dn['e'] += rate_ion
        
        dE['H'] -= rate_ion * TH
        dE['Hi'] += rate_ion * TH
        dE['e'] -= rate_ion * 13.6  # Ionization energy loss
        
        # 3. Three-body recombination: e + e + H+ -> e + H
        # Note: 3-body rate is cm⁶/s, need CM3_TO_M3² = 1e-12
        k_3body = rates.H_3body_recombination(Te, rates.H_ionization(Te)) * ndamp(nHi_cgs) * CM3_TO_M3 * CM3_TO_M3
        rate_3body = safe_rate(k_3body, ne, ne, nHi)
        
        dn['Hi'] -= rate_3body
        dn['H'] += rate_3body
        dn['e'] -= rate_3body
        
        dE['Hi'] -= rate_3body * THi
        dE['H'] += rate_3body * THi
        dE['e'] += rate_3body * (13.6 - Te / 3.0)  # Energy gain
        
        # 4. Radiative recombination: e + H+ -> H + photon
        k_rad = rates.H_radiative_recombination(Te) * ndamp(nHi_cgs) * CM3_TO_M3
        rate_rad = safe_rate(k_rad, ne, nHi)
        
        dn['Hi'] -= rate_rad
        dn['e'] -= rate_rad
        dn['H'] += rate_rad
        
        dE['Hi'] -= rate_rad * THi
        dE['H'] += rate_rad * THi
        # Electron energy loss with HYDHEL correction
        dE['e'] -= rate_rad * Te * 0.667 * 1.5  # Approximate
    
    # -------------------------------------------------------------------
    # Helium Reactions (bHe)
    # Rate coefficients from rates.py are in cm³/s (with 2D ADAS interpolation)
    # Need to pass ne in CGS for 2D interpolation
    # ndamp() takes CGS densities
    # -------------------------------------------------------------------
    if include_He and nHeI is not None:
        # Convert He densities to CGS
        nHeI_cgs = nHeI * M3_TO_CM3
        nHeII_cgs = nHeII * M3_TO_CM3
        nHeIII_cgs = nHeIII * M3_TO_CM3
        
        # HeI ionization: e + He -> e + He+ + e
        # Pass ne_cgs for 2D ADAS interpolation
        k_HeI_ion = rates.HeI_ionization(Te, ne_cgs) * ndamp(nHeI_cgs) * CM3_TO_M3
        rate_HeI_ion = safe_rate(k_HeI_ion, ne, nHeI)
        
        dn['HeI'] -= rate_HeI_ion
        dn['HeII'] += rate_HeI_ion
        dn['e'] += rate_HeI_ion
        
        dE['HeI'] -= rate_HeI_ion * THeI
        dE['HeII'] += rate_HeI_ion * THeI
        
        # HeII ionization: e + He+ -> e + He++ + e
        k_HeII_ion = rates.HeII_ionization(Te, ne_cgs) * ndamp(nHeII_cgs) * CM3_TO_M3
        rate_HeII_ion = safe_rate(k_HeII_ion, ne, nHeII)
        
        dn['HeII'] -= rate_HeII_ion
        dn['HeIII'] += rate_HeII_ion
        dn['e'] += rate_HeII_ion
        
        dE['HeII'] -= rate_HeII_ion * THeII
        dE['HeIII'] += rate_HeII_ion * THeII
        
        # Helium cooling (excitation + ionization energy loss)
        # Cooling rates L are in [eV*m³/s], so L * n_He * n_e gives [eV/m³/s]
        if use_ADAS:
            L_HeI, L_HeII, L_HeIII = rates.He_cooling_rate(Te)
            # ADAS rates in eV*m³/s
            cooling = safe_rate(L_HeI, nHeI, ne) + safe_rate(L_HeII, nHeII, ne) + safe_rate(L_HeIII, nHeIII, ne)
            dE['e'] -= (2.0/3.0) * cooling
        else:
            L_HeI, L_HeII, L_HeIII = rates.He_cooling_rate_IAEA(Te)
            # IAEA rates in eV*cm³/s, convert to eV*m³/s: multiply by CM3_TO_M3
            cooling = safe_rate(L_HeI * CM3_TO_M3, nHeI, ne) + safe_rate(L_HeII * CM3_TO_M3, nHeII, ne) + safe_rate(L_HeIII * CM3_TO_M3, nHeIII, ne)
            dE['e'] -= (2.0/3.0) * cooling
        
        # HeII recombination: e + He+ -> He + photon
        k_HeII_rec = rates.HeII_recombination(Te, ne_cgs) * ndamp(nHeII_cgs) * CM3_TO_M3
        rate_HeII_rec = safe_rate(k_HeII_rec, ne, nHeII)
        
        dn['HeI'] += rate_HeII_rec
        dn['HeII'] -= rate_HeII_rec
        dn['e'] -= rate_HeII_rec
        
        dE['HeI'] += rate_HeII_rec * THeII
        dE['HeII'] -= rate_HeII_rec * THeII
        dE['e'] -= rate_HeII_rec * Te * 0.667 * 1.5
        
        # HeIII recombination: e + He++ -> He+ + photon
        k_HeIII_rec = rates.HeIII_recombination(Te, ne_cgs) * ndamp(nHeIII_cgs) * CM3_TO_M3
        rate_HeIII_rec = safe_rate(k_HeIII_rec, ne, nHeIII)
        
        dn['HeII'] += rate_HeIII_rec
        dn['HeIII'] -= rate_HeIII_rec
        dn['e'] -= rate_HeIII_rec
        
        dE['HeII'] += rate_HeIII_rec * THeIII
        dE['HeIII'] -= rate_HeIII_rec * THeIII
        dE['e'] -= rate_HeIII_rec * Te * 0.667 * 1.5
    
    # -------------------------------------------------------------------
    # Charge Exchange Reactions (bcx)
    # Rate coefficients from rates.py are in cm³/s, convert to m³/s
    # -------------------------------------------------------------------
    if include_cx:
        # H+ + H -> H + H+ (energy exchange only, no particle change)
        k_HiH_cx = rates.HiH_charge_exchange(THi, TH) * ndamp(nHi) * CM3_TO_M3
        rate_HiH_cx = safe_rate(k_HiH_cx, nHi, nH)
        
        dE['Hi'] += rate_HiH_cx * (TH - THi)
        dE['H'] += rate_HiH_cx * (THi - TH)
        
        # Helium charge exchange reactions
        if include_He and nHeI is not None:
            # He+ + H -> He + H+
            k_HeIIH_cx = rates.HeIIH_charge_exchange(THeII, TH) * ndamp(nHeII) * CM3_TO_M3
            rate_HeIIH_cx = safe_rate(k_HeIIH_cx, nHeII, nH)
            
            dn['HeII'] -= rate_HeIIH_cx
            dn['HeI'] += rate_HeIIH_cx
            dn['H'] -= rate_HeIIH_cx
            dn['Hi'] += rate_HeIIH_cx
            
            dE['HeII'] -= rate_HeIIH_cx * THeII
            dE['HeI'] += rate_HeIIH_cx * THeII
            dE['H'] -= rate_HeIIH_cx * TH
            dE['Hi'] += rate_HeIIH_cx * TH
            
            # He++ + H -> He+ + H+
            k_HeIIIH_cx = rates.HeIIIH_charge_exchange(THeIII, TH) * ndamp(nHeIII) * CM3_TO_M3
            rate_HeIIIH_cx = safe_rate(k_HeIIIH_cx, nHeIII, nH)
            
            dn['HeIII'] -= rate_HeIIIH_cx
            dn['HeII'] += rate_HeIIIH_cx
            dn['H'] -= rate_HeIIIH_cx
            dn['Hi'] += rate_HeIIIH_cx
            
            dE['HeIII'] -= rate_HeIIIH_cx * THeIII
            dE['HeII'] += rate_HeIIIH_cx * THeIII
            dE['H'] -= rate_HeIIIH_cx * TH
            dE['Hi'] += rate_HeIIIH_cx * TH
    
    return dn, dE


def compute_sources_from_state(state, params: dict) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """
    Compute collision sources from a PlasmaState object.
    
    Parameters
    ----------
    state : PlasmaState
        Current plasma state.
    params : dict
        Simulation parameters with physics flags.
        
    Returns
    -------
    dn : Dict[str, np.ndarray]
        Density source terms [m^-3/s].
    dE : Dict[str, np.ndarray]
        Energy source terms [eV·m^-3/s].
    """
    # Extract density and temperature arrays
    ne = state.electrons.n.x.array.copy()
    Te = state.electrons.T.copy()
    
    nH = state.species['H'].n.x.array.copy() if 'H' in state.species else np.zeros_like(ne)
    TH = state.species['H'].T.copy() if 'H' in state.species else np.ones_like(ne)
    
    nHi = state.species['Hi'].n.x.array.copy() if 'Hi' in state.species else np.zeros_like(ne)
    THi = state.species['Hi'].T.copy() if 'Hi' in state.species else Te.copy()
    
    # Helium species (optional)
    nHeI = state.species['HeI'].n.x.array.copy() if 'HeI' in state.species else None
    THeI = state.species['HeI'].T.copy() if 'HeI' in state.species else None
    
    nHeII = state.species['HeII'].n.x.array.copy() if 'HeII' in state.species else None
    THeII = state.species['HeII'].T.copy() if 'HeII' in state.species else None
    
    nHeIII = state.species['HeIII'].n.x.array.copy() if 'HeIII' in state.species else None
    THeIII = state.species['HeIII'].T.copy() if 'HeIII' in state.species else None
    
    # Get physics flags from params
    include_H = params.get('bH', True)
    include_He = params.get('bHe', 'HeI' in state.species)
    include_cx = params.get('bcx', True)
    use_ADAS = params.get('bADAS', True)
    
    return compute_collision_sources(
        ne, Te, nH, TH, nHi, THi,
        nHeI, THeI, nHeII, THeII, nHeIII, THeIII,
        use_ADAS=use_ADAS,
        include_H=include_H,
        include_He=include_He,
        include_cx=include_cx
    )
