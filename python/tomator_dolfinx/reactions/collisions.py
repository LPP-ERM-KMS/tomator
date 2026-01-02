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

# Temperature clamps matching C++ Tomator1D (reactionrates.cpp)
# RRH2 uses Te_max = 100 eV, RRH3 uses Te_max = 1000 eV, RR uses T_max = 2e4 eV
# Heavy particle temps should be bounded similar to RRH3
T_MIN = 0.026  # Room temperature ~0.026 eV
T_MAX_HEAVY = 1000.0  # Max for ions/neutrals (matching C++ RRH3)
T_MAX_ELECTRON = 2e4  # Max for electrons (matching C++ RR)


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
    nH2: np.ndarray = None,
    TH2: np.ndarray = None,
    nH2i: np.ndarray = None,
    TH2i: np.ndarray = None,
    nH3i: np.ndarray = None,
    TH3i: np.ndarray = None,
    nHeI: np.ndarray = None,
    THeI: np.ndarray = None,
    nHeII: np.ndarray = None,
    THeII: np.ndarray = None,
    nHeIII: np.ndarray = None,
    THeIII: np.ndarray = None,
    use_ADAS: bool = True,
    include_H: bool = True,
    include_H2: bool = True,
    include_He: bool = True,
    include_ion: bool = True,
    include_cx: bool = True,
    include_elastic: bool = True,
    include_coulomb: bool = True
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """
    Compute all collision source terms for density, energy, and collision frequency.
    
    This is the main function that computes dn/dt, dE/dt, and nu for each species
    from all relevant reactions. The collision frequencies are used for computing
    physics-based transport coefficients (diffusion).
    
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
    nH2 : np.ndarray, optional
        Molecular hydrogen density [m^-3].
    TH2 : np.ndarray, optional
        Molecular hydrogen temperature [eV].
    nH2i : np.ndarray, optional
        H2+ density [m^-3].
    TH2i : np.ndarray, optional
        H2+ temperature [eV].
    nH3i : np.ndarray, optional
        H3+ density [m^-3].
    TH3i : np.ndarray, optional
        H3+ temperature [eV].
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
        Include atomic hydrogen reactions.
    include_H2 : bool
        Include molecular hydrogen (H2, H2+, H3+) reactions.
    include_He : bool
        Include helium reactions.
    include_ion : bool
        Include ion-neutral reactions (proton impact excitation, ionization, etc.).
    include_cx : bool
        Include charge exchange reactions.
    include_elastic : bool
        Include elastic collision energy exchange.
    include_coulomb : bool
        Include Coulomb collision energy exchange between charged particles.
        
    Returns
    -------
    dn : Dict[str, np.ndarray]
        Density source terms [m^-3/s].
    dE : Dict[str, np.ndarray]
        Energy source terms [eV·m^-3/s].
    nu : Dict[str, np.ndarray]
        Collision frequencies [1/s] for each species.
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
    # Collision frequency arrays [1/s]
    nu = {
        'e': np.zeros(nmesh),
        'H': np.zeros(nmesh),
        'Hi': np.zeros(nmesh),
    }
    
    # Add H2 species if provided
    if nH2 is not None:
        dn['H2'] = np.zeros(nmesh)
        dn['H2i'] = np.zeros(nmesh)
        dn['H3i'] = np.zeros(nmesh)
        dE['H2'] = np.zeros(nmesh)
        dE['H2i'] = np.zeros(nmesh)
        dE['H3i'] = np.zeros(nmesh)
        nu['H2'] = np.zeros(nmesh)
        nu['H2i'] = np.zeros(nmesh)
        nu['H3i'] = np.zeros(nmesh)
    
    # Add He species if provided
    if nHeI is not None:
        dn['HeI'] = np.zeros(nmesh)
        dn['HeII'] = np.zeros(nmesh)
        dn['HeIII'] = np.zeros(nmesh)
        dE['HeI'] = np.zeros(nmesh)
        dE['HeII'] = np.zeros(nmesh)
        dE['HeIII'] = np.zeros(nmesh)
        nu['HeI'] = np.zeros(nmesh)
        nu['HeII'] = np.zeros(nmesh)
        nu['HeIII'] = np.zeros(nmesh)
    
    rates = ReactionRates()
    
    # -------------------------------------------------------------------
    # Clamp temperatures to prevent unphysical values (matching C++ Tomator1D)
    # Heavy particle temps bounded to 1000 eV (RRH3 limit)
    # Electron temp bounded to 2e4 eV (HYDHEL limit)
    # -------------------------------------------------------------------
    Te = np.clip(Te, T_MIN, T_MAX_ELECTRON)
    TH = np.clip(TH, T_MIN, T_MAX_HEAVY)
    THi = np.clip(THi, T_MIN, T_MAX_HEAVY)
    if TH2 is not None:
        TH2 = np.clip(TH2, T_MIN, T_MAX_HEAVY)
    if TH2i is not None:
        TH2i = np.clip(TH2i, T_MIN, T_MAX_HEAVY)
    if TH3i is not None:
        TH3i = np.clip(TH3i, T_MIN, T_MAX_HEAVY)
    if THeI is not None:
        THeI = np.clip(THeI, T_MIN, T_MAX_HEAVY)
    if THeII is not None:
        THeII = np.clip(THeII, T_MIN, T_MAX_HEAVY)
    if THeIII is not None:
        THeIII = np.clip(THeIII, T_MIN, T_MAX_HEAVY)
    
    # Convert SI densities to CGS for rate calculations (m^-3 -> cm^-3)
    ne_cgs = ne * M3_TO_CM3
    nH_cgs = nH * M3_TO_CM3
    nHi_cgs = nHi * M3_TO_CM3
    
    # -------------------------------------------------------------------
    # Electron-Hydrogen Reactions (bH)
    # Rate coefficients from rates.py are in cm³/s, convert to m³/s
    # ndamp() takes CGS densities
    # nu_X += k * n_partner [1/s] (like C++ Tomator1D collisions.cpp)
    # -------------------------------------------------------------------
    if include_H:
        # 1. H excitation: e + H -> e + H* (energy loss only)
        k_exc = rates.H_excitation(Te) * CM3_TO_M3
        rate_exc = safe_rate(k_exc, ne, nH)
        dE['e'] -= rate_exc * 10.2  # Excitation energy loss
        # Collision frequency: H collides with electrons
        nu['e'] += k_exc * nH
        nu['H'] += k_exc * ne
        
        # 2. H ionization: e + H -> e + H+ + e
        k_ion = rates.H_ionization(Te) * ndamp(nH_cgs) * CM3_TO_M3
        rate_ion = safe_rate(k_ion, ne, nH)
        
        dn['H'] -= rate_ion
        dn['Hi'] += rate_ion
        dn['e'] += rate_ion
        
        dE['H'] -= rate_ion * TH
        dE['Hi'] += rate_ion * TH
        dE['e'] -= rate_ion * 13.6  # Ionization energy loss
        # Collision frequency
        nu['e'] += k_ion * nH
        nu['H'] += k_ion * ne
        
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
        # 3-body nu: k_3body * ne * nHi for electrons
        nu['e'] += k_3body * ne * nHi
        nu['Hi'] += k_3body * ne * ne
        
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
        # Collision frequency
        nu['e'] += k_rad * nHi
        nu['Hi'] += k_rad * ne
    
    # -------------------------------------------------------------------
    # Molecular Hydrogen Reactions (bH2)
    # Rate coefficients from rates.py using Dirk W.'s 2D tables
    # Energy thresholds from collisions.cpp
    # -------------------------------------------------------------------
    if include_H2 and nH2 is not None:
        # Convert H2 densities to CGS
        nH2_cgs = nH2 * M3_TO_CM3
        nH2i_cgs = nH2i * M3_TO_CM3 if nH2i is not None else np.zeros_like(nH2_cgs)
        nH3i_cgs = nH3i * M3_TO_CM3 if nH3i is not None else np.zeros_like(nH2_cgs)
        
        # Ensure H2i and H3i arrays exist
        if nH2i is None:
            nH2i = np.zeros(nmesh)
            TH2i = Te.copy()
        if nH3i is None:
            nH3i = np.zeros(nmesh)
            TH3i = Te.copy()
        
        # Energy thresholds from collisions.cpp
        E_H2_diss = 4.52       # H2 dissociation energy [eV]
        E_H2_ion = 15.43       # H2 ionization energy [eV]
        E_H2i_diss_rec = 10.5  # H2+ dissociative recombination energy release [eV]
        
        # ----- e + H2 elastic collision (energy exchange) -----
        # From collisions.cpp: uses RRH2(ELAS, ne, Te) with Langevin factor
        # L = 4*me*mH2/(me+mH2)^2 ≈ 4*me*2mi/(me+2mi)^2 ≈ 1.09e-3
        me_over_mH2 = 1.0 / 3672.0  # me/mH2 ratio
        L_eH2 = 4.0 * me_over_mH2 / (1.0 + me_over_mH2)**2  # ~1.09e-3
        k_eH2_elas = rates.H2_elastic(Te) * CM3_TO_M3
        rate_eH2_elas = safe_rate(k_eH2_elas, ne, nH2)
        dE['e'] += rate_eH2_elas * L_eH2 * (TH2 - Te)
        dE['H2'] += rate_eH2_elas * L_eH2 * (Te - TH2)
        # Collision frequencies
        nu['e'] += k_eH2_elas * nH2
        nu['H2'] += k_eH2_elas * ne
        
        # ----- Reaction 2.2.1a: e + H2(v=0) -> e + H2(v=1) vibrational excitation -----
        # Energy loss: 0.5 eV (no density change)
        mask_exc = Te > 0.1
        k_H2_exc_a = rates.H2_excitation_vib_a(Te) * CM3_TO_M3
        rate_H2_exc_a = safe_rate(k_H2_exc_a, ne, nH2)
        rate_H2_exc_a = np.where(mask_exc, rate_H2_exc_a, 0.0)
        dE['e'] -= rate_H2_exc_a * 0.5
        # nu for vibrational excitation
        nu['e'] += np.where(mask_exc, k_H2_exc_a * nH2, 0.0)
        nu['H2'] += np.where(mask_exc, k_H2_exc_a * ne, 0.0)
        
        # ----- Reaction 2.2.1b: e + H2(v=0) -> e + H2(v=2) vibrational excitation -----
        # Energy loss: 1.0 eV (no density change)
        k_H2_exc_b = rates.H2_excitation_vib_b(Te) * CM3_TO_M3
        rate_H2_exc_b = safe_rate(k_H2_exc_b, ne, nH2)
        rate_H2_exc_b = np.where(mask_exc, rate_H2_exc_b, 0.0)
        dE['e'] -= rate_H2_exc_b * 1.0
        nu['e'] += np.where(mask_exc, k_H2_exc_b * nH2, 0.0)
        nu['H2'] += np.where(mask_exc, k_H2_exc_b * ne, 0.0)
        
        # ----- Reaction 2.2.2: e + H2(X) -> e + H2(B) electronic excitation -----
        # Energy loss: 12.1 eV (no density change)
        k_H2_exc_B = rates.H2_excitation_elec_B(Te) * CM3_TO_M3
        rate_H2_exc_B = safe_rate(k_H2_exc_B, ne, nH2)
        dE['e'] -= rate_H2_exc_B * 12.1
        nu['e'] += k_H2_exc_B * nH2
        nu['H2'] += k_H2_exc_B * ne
        
        # ----- Reaction 2.2.3: e + H2(X) -> e + H2(C) electronic excitation -----
        # Energy loss: 12.4 eV (no density change)
        k_H2_exc_C = rates.H2_excitation_elec_C(Te) * CM3_TO_M3
        rate_H2_exc_C = safe_rate(k_H2_exc_C, ne, nH2)
        dE['e'] -= rate_H2_exc_C * 12.4
        nu['e'] += k_H2_exc_C * nH2
        nu['H2'] += k_H2_exc_C * ne
        
        # ----- Reaction 2.2.4: e + H2(X) -> e + H2(E,F) electronic excitation -----
        # Energy loss: 12.7 eV (no density change)
        k_H2_exc_EF = rates.H2_excitation_elec_EF(Te) * CM3_TO_M3
        rate_H2_exc_EF = safe_rate(k_H2_exc_EF, ne, nH2)
        dE['e'] -= rate_H2_exc_EF * 12.7
        nu['e'] += k_H2_exc_EF * nH2
        nu['H2'] += k_H2_exc_EF * ne
        
        # ----- Reaction 2.2.5-2.2.8: e + H2 -> e + H + H (dissociation) -----
        # Note: This is the actual dissociation using RRH2(DISS,...) 2D table
        k_H2_diss = rates.H2_dissociation(Te, ne_cgs) * ndamp(nH2_cgs) * CM3_TO_M3
        rate_H2_diss = safe_rate(k_H2_diss, ne, nH2)
        
        dn['H2'] -= rate_H2_diss
        dn['H'] += 2.0 * rate_H2_diss  # Creates 2 H atoms
        
        # From C++: dEe -= 10.5, dEH2 -= TH2, dEH += (TH2 + 2*3.0)
        dE['e'] -= rate_H2_diss * 10.5
        dE['H2'] -= rate_H2_diss * TH2
        dE['H'] += rate_H2_diss * (TH2 + 6.0)  # TH2 + 2*3.0 eV kinetic energy
        # Collision frequencies
        nu['e'] += k_H2_diss * nH2
        nu['H2'] += k_H2_diss * ne
        
        # ----- Reaction 2.2.9: e + H2 -> 2e + H2+ (ionization) -----
        k_H2_ion = rates.H2_ionization(Te, ne_cgs) * ndamp(nH2_cgs) * CM3_TO_M3
        rate_H2_ion = safe_rate(k_H2_ion, ne, nH2)
        
        dn['H2'] -= rate_H2_ion
        dn['H2i'] += rate_H2_ion
        dn['e'] += rate_H2_ion
        
        dE['H2'] -= rate_H2_ion * TH2
        dE['H2i'] += rate_H2_ion * TH2
        dE['e'] -= rate_H2_ion * E_H2_ion  # Electron energy loss
        # Collision frequencies
        nu['e'] += k_H2_ion * nH2
        nu['H2'] += k_H2_ion * ne
        
        # ----- RECO: e + H2+ -> H2 (recombination to ground state H2) -----
        # From C++ collisions.cpp lines 306-317: uses RRH2(RECO) and produces H2
        k_H2i_rec = rates.H2i_recombination(Te) * ndamp(nH2i_cgs) * CM3_TO_M3
        rate_H2i_rec = safe_rate(k_H2i_rec, ne, nH2i)
        
        dn['H2i'] -= rate_H2i_rec
        dn['e'] -= rate_H2i_rec
        dn['H2'] += rate_H2i_rec  # Produces H2 (not 2H!)
        
        # From C++: dEe -= knn * Te * 0.667 * 0.89
        dE['e'] -= rate_H2i_rec * Te * 0.667 * 0.89  # HYDHEL correction
        dE['H2i'] -= rate_H2i_rec * TH2i
        dE['H2'] += rate_H2i_rec * TH2i
        # Collision frequencies
        nu['e'] += k_H2i_rec * nH2i
        nu['H2i'] += k_H2i_rec * ne
        
        # ----- Reaction 2.2.14: e + H2+ -> H + H (dissociative recombination) -----
        # From C++ collisions.cpp lines 373-387: uses RR(REAC2214) and produces 2H atoms
        k_H2i_diss_rec = rates.H2i_dissociative_recombination(Te) * ndamp(nH2i_cgs) * CM3_TO_M3
        rate_H2i_diss_rec = safe_rate(k_H2i_diss_rec, ne, nH2i)
        
        dn['H2i'] -= rate_H2i_diss_rec
        dn['e'] -= rate_H2i_diss_rec
        dn['H'] += 2.0 * rate_H2i_diss_rec  # Creates 2 H atoms
        
        # From C++: dEe -= knn * Te * 0.667 * (1.5 + d(ln k)/d(ln Te))
        # Approximate the derivative term as ~0.5 for simplicity
        dE['e'] -= rate_H2i_diss_rec * Te * 0.667 * 2.0  # ~(1.5 + 0.5)
        dE['H2i'] -= rate_H2i_diss_rec * TH2i
        dE['H'] += rate_H2i_diss_rec * TH2i  # Energy goes to H atoms
        # Collision frequencies
        nu['e'] += k_H2i_diss_rec * nH2i
        nu['H2i'] += k_H2i_diss_rec * ne
        
        # ----- Reaction 2.2.10: e + H2 -> e + H+ + H + e (dissociative ionization of H2) -----
        # From C++: dnH2 -= knn, dne += knn, dnHi += knn, dnH += knn
        # Energy: dEe -= 18.0, dEH2 -= TH2, dEHi += TH2/2 + 0.1, dEH += TH2/2 + 0.1
        k_H2_diss_ion = rates.H2i_dissociative_ionization(Te) * ndamp(nH2_cgs) * CM3_TO_M3
        rate_H2_diss_ion = safe_rate(k_H2_diss_ion, ne, nH2)
        
        dn['H2'] -= rate_H2_diss_ion
        dn['e'] += rate_H2_diss_ion  # Net +1 electron
        dn['Hi'] += rate_H2_diss_ion
        dn['H'] += rate_H2_diss_ion
        
        dE['H2'] -= rate_H2_diss_ion * TH2
        dE['Hi'] += rate_H2_diss_ion * (TH2 / 2.0 + 0.1)
        dE['H'] += rate_H2_diss_ion * (TH2 / 2.0 + 0.1)
        dE['e'] -= rate_H2_diss_ion * 18.0  # Energy cost from C++
        # Collision frequencies
        nu['e'] += k_H2_diss_ion * nH2
        nu['H2'] += k_H2_diss_ion * ne
        
        # ----- Reaction 2.2.15a: e + H3+ -> H2 + H -----
        # From C++ collisions.cpp lines 404-416: uses REAC2215 for both H3+ recombination channels
        k_H3i_rec_H2H = rates.H3i_recombination_3H(Te) * ndamp(nH3i_cgs) * CM3_TO_M3  # Same rate as 3H channel
        rate_H3i_rec_H2H = safe_rate(k_H3i_rec_H2H, ne, nH3i)
        
        dn['H3i'] -= rate_H3i_rec_H2H
        dn['e'] -= rate_H3i_rec_H2H
        dn['H2'] += rate_H3i_rec_H2H
        dn['H'] += rate_H3i_rec_H2H
        
        # From C++: dEe -= knn * Te * 0.667 * (~2), dEH2 += (TH3i + Te/3) * 2/3, dEH += (TH3i + Te/3) * 1/3
        dE['e'] -= rate_H3i_rec_H2H * Te * 0.667 * 2.0  # HYDHEL correction
        dE['H3i'] -= rate_H3i_rec_H2H * TH3i
        dE['H2'] += rate_H3i_rec_H2H * (TH3i + Te / 3.0) * 2.0 / 3.0
        dE['H'] += rate_H3i_rec_H2H * (TH3i + Te / 3.0) * 1.0 / 3.0
        # Collision frequencies
        nu['e'] += k_H3i_rec_H2H * nH3i
        nu['H3i'] += k_H3i_rec_H2H * ne
        
        # ----- Reaction 2.2.15b: e + H3+ -> H + H + H -----
        # From C++ collisions.cpp lines 389-401
        k_H3i_rec_3H = rates.H3i_recombination_3H(Te) * ndamp(nH3i_cgs) * CM3_TO_M3
        rate_H3i_rec_3H = safe_rate(k_H3i_rec_3H, ne, nH3i)
        
        dn['H3i'] -= rate_H3i_rec_3H
        dn['e'] -= rate_H3i_rec_3H
        dn['H'] += 3.0 * rate_H3i_rec_3H
        
        # From C++: dEe -= knn * Te * 0.667 * (~2), dEH += (TH3i + Te/3)
        dE['e'] -= rate_H3i_rec_3H * Te * 0.667 * 2.0  # HYDHEL correction
        dE['H3i'] -= rate_H3i_rec_3H * TH3i
        dE['H'] += rate_H3i_rec_3H * (TH3i + Te / 3.0)
        # Collision frequencies
        nu['e'] += k_H3i_rec_3H * nH3i
        nu['H3i'] += k_H3i_rec_3H * ne
        
        # ----- Reaction 2.2.12: e + H2+ -> e + H+ + H (dissociation) -----
        # From collisions.cpp: dEe -= 10.5, dEHi += (TH2i/2 + 4.3), dEH += (TH2i/2 + 4.3)
        k_H2i_diss = rates.H2i_dissociation(Te) * ndamp(nH2i_cgs) * CM3_TO_M3
        rate_H2i_diss = safe_rate(k_H2i_diss, ne, nH2i)
        
        dn['H2i'] -= rate_H2i_diss
        dn['Hi'] += rate_H2i_diss
        dn['H'] += rate_H2i_diss
        
        dE['e'] -= rate_H2i_diss * 10.5
        dE['H2i'] -= rate_H2i_diss * TH2i
        dE['Hi'] += rate_H2i_diss * (TH2i / 2.0 + 4.3)
        dE['H'] += rate_H2i_diss * (TH2i / 2.0 + 4.3)
        # Collision frequencies
        nu['e'] += k_H2i_diss * nH2i
        nu['H2i'] += k_H2i_diss * ne
        
        # ----- Reaction 2.2.13: e + H2+ -> e + H+ + H* (dissociation with excitation) -----
        # From collisions.cpp: dEe -= 17.5, dEHi += (TH2i/2 + 1.5), dEH += (TH2i/2 + 1.5)
        k_H2i_diss_exc = rates.H2i_dissociation_excitation(Te) * ndamp(nH2i_cgs) * CM3_TO_M3
        rate_H2i_diss_exc = safe_rate(k_H2i_diss_exc, ne, nH2i)
        
        dn['H2i'] -= rate_H2i_diss_exc
        dn['Hi'] += rate_H2i_diss_exc
        dn['H'] += rate_H2i_diss_exc
        
        dE['e'] -= rate_H2i_diss_exc * 17.5
        dE['H2i'] -= rate_H2i_diss_exc * TH2i
        dE['Hi'] += rate_H2i_diss_exc * (TH2i / 2.0 + 1.5)
        dE['H'] += rate_H2i_diss_exc * (TH2i / 2.0 + 1.5)
        # Collision frequencies
        nu['e'] += k_H2i_diss_exc * nH2i
        nu['H2i'] += k_H2i_diss_exc * ne
        
        # ----- Reaction 2.2.16: e + H3+ -> e + H+ + 2H (dissociation) -----
        # From collisions.cpp: dEe -= 14.0, dEHi += (TH3i/3 + 4.33), dEH += 2*(TH3i/3 + 4.33)
        k_H3i_diss = rates.H3i_dissociation(Te) * ndamp(nH3i_cgs) * CM3_TO_M3
        rate_H3i_diss = safe_rate(k_H3i_diss, ne, nH3i)
        
        dn['H3i'] -= rate_H3i_diss
        dn['Hi'] += rate_H3i_diss
        dn['H'] += 2.0 * rate_H3i_diss
        
        dE['e'] -= rate_H3i_diss * 14.0
        dE['H3i'] -= rate_H3i_diss * TH3i
        dE['Hi'] += rate_H3i_diss * (TH3i / 3.0 + 4.33)
        dE['H'] += rate_H3i_diss * 2.0 * (TH3i / 3.0 + 4.33)
        # Collision frequencies
        nu['e'] += k_H3i_diss * nH3i
        nu['H3i'] += k_H3i_diss * ne
        
        # ----- Reaction 3.2.3: H+ + H2 -> H2+ + H (charge exchange) -----
        k_Hi_H2_cx = rates.Hi_H2_charge_exchange(THi, TH2) * ndamp(nHi_cgs) * CM3_TO_M3
        rate_Hi_H2_cx = safe_rate(k_Hi_H2_cx, nHi, nH2)
        
        dn['Hi'] -= rate_Hi_H2_cx
        dn['H2'] -= rate_Hi_H2_cx
        dn['H2i'] += rate_Hi_H2_cx
        dn['H'] += rate_Hi_H2_cx
        
        dE['Hi'] -= rate_Hi_H2_cx * THi
        dE['H2'] -= rate_Hi_H2_cx * TH2
        dE['H2i'] += rate_Hi_H2_cx * THi  # H2i gets H+ energy
        dE['H'] += rate_Hi_H2_cx * TH2     # H gets H2 energy
        # Collision frequencies
        nu['Hi'] += k_Hi_H2_cx * nH2
        nu['H2'] += k_Hi_H2_cx * nHi
        
        # ----- Reaction 3.2.8: H+ + H2 -> H+ + H + H (dissociation) -----
        k_Hi_H2_diss = rates.Hi_H2_dissociation(THi, TH2) * ndamp(nH2_cgs) * CM3_TO_M3
        rate_Hi_H2_diss = safe_rate(k_Hi_H2_diss, nHi, nH2)
        
        dn['H2'] -= rate_Hi_H2_diss
        dn['H'] += 2.0 * rate_Hi_H2_diss
        
        dE['H2'] -= rate_Hi_H2_diss * TH2
        dE['H'] += rate_Hi_H2_diss * TH2
        dE['Hi'] -= rate_Hi_H2_diss * E_H2_diss  # Energy comes from H+
        # Collision frequencies
        nu['Hi'] += k_Hi_H2_diss * nH2
        nu['H2'] += k_Hi_H2_diss * nHi
        
        # ----- Reaction 4.2.5: H2+ + H2 -> H3+ + H -----
        k_H2i_H2_H3i = rates.H2i_H2_to_H3i(TH2i, TH2) * ndamp(nH2i_cgs) * CM3_TO_M3
        rate_H2i_H2_H3i = safe_rate(k_H2i_H2_H3i, nH2i, nH2)
        
        dn['H2i'] -= rate_H2i_H2_H3i
        dn['H2'] -= rate_H2i_H2_H3i
        dn['H3i'] += rate_H2i_H2_H3i
        dn['H'] += rate_H2i_H2_H3i
        
        dE['H2i'] -= rate_H2i_H2_H3i * TH2i
        dE['H2'] -= rate_H2i_H2_H3i * TH2
        dE['H3i'] += rate_H2i_H2_H3i * (TH2i + TH2) * 0.75  # 3/4 to H3+
        dE['H'] += rate_H2i_H2_H3i * (TH2i + TH2) * 0.25    # 1/4 to H
        # Collision frequencies
        nu['H2i'] += k_H2i_H2_H3i * nH2
        nu['H2'] += k_H2i_H2_H3i * nH2i
        
        # ----- Reaction 4.1.9: H2 + H -> H + H + H (molecular dissociation by H impact) -----
        k_H2_H_diss = rates.H2_H_dissociation(TH2, TH) * ndamp(nH2_cgs) * CM3_TO_M3
        rate_H2_H_diss = safe_rate(k_H2_H_diss, nH2, nH)
        
        dn['H2'] -= rate_H2_H_diss
        dn['H'] += 2.0 * rate_H2_H_diss  # Net +2 H atoms (H2->2H, H->H)
        
        dE['H2'] -= rate_H2_H_diss * TH2
        dE['H'] += rate_H2_H_diss * TH2  # Energy transfer
        dE['H'] -= rate_H2_H_diss * E_H2_diss  # Energy cost
        # Collision frequencies
        nu['H2'] += k_H2_H_diss * nH
        nu['H'] += k_H2_H_diss * nH2
    
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
        # Collision frequencies
        nu['e'] += k_HeI_ion * nHeI
        nu['HeI'] += k_HeI_ion * ne
        
        # HeII ionization: e + He+ -> e + He++ + e
        k_HeII_ion = rates.HeII_ionization(Te, ne_cgs) * ndamp(nHeII_cgs) * CM3_TO_M3
        rate_HeII_ion = safe_rate(k_HeII_ion, ne, nHeII)
        
        dn['HeII'] -= rate_HeII_ion
        dn['HeIII'] += rate_HeII_ion
        dn['e'] += rate_HeII_ion
        
        dE['HeII'] -= rate_HeII_ion * THeII
        dE['HeIII'] += rate_HeII_ion * THeII
        # Collision frequencies
        nu['e'] += k_HeII_ion * nHeII
        nu['HeII'] += k_HeII_ion * ne
        
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
        # Collision frequencies
        nu['e'] += k_HeII_rec * nHeII
        nu['HeII'] += k_HeII_rec * ne
        
        # HeIII recombination: e + He++ -> He+ + photon
        k_HeIII_rec = rates.HeIII_recombination(Te, ne_cgs) * ndamp(nHeIII_cgs) * CM3_TO_M3
        rate_HeIII_rec = safe_rate(k_HeIII_rec, ne, nHeIII)
        
        dn['HeII'] += rate_HeIII_rec
        dn['HeIII'] -= rate_HeIII_rec
        dn['e'] -= rate_HeIII_rec
        
        dE['HeII'] += rate_HeIII_rec * THeIII
        dE['HeIII'] -= rate_HeIII_rec * THeIII
        dE['e'] -= rate_HeIII_rec * Te * 0.667 * 1.5
        # Collision frequencies
        nu['e'] += k_HeIII_rec * nHeIII
        nu['HeIII'] += k_HeIII_rec * ne
    
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
        # Collision frequencies
        nu['Hi'] += k_HiH_cx * nH
        nu['H'] += k_HiH_cx * nHi
        
        # H2+ + H2 -> H2 + H2+ (symmetric charge exchange, energy exchange only)
        if include_H2:
            nH2i_cgs = nH2i * M3_TO_CM3 if nH2i is not None else 0.0
            k_H2iH2_cx = rates.H2iH2_charge_exchange(TH2i, TH2) * ndamp(nH2i_cgs) * CM3_TO_M3
            rate_H2iH2_cx = safe_rate(k_H2iH2_cx, nH2i, nH2)
            
            dE['H2i'] += rate_H2iH2_cx * (TH2 - TH2i)
            dE['H2'] += rate_H2iH2_cx * (TH2i - TH2)
            # Collision frequencies
            nu['H2i'] += k_H2iH2_cx * nH2
            nu['H2'] += k_H2iH2_cx * nH2i
        
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
            # Collision frequencies
            nu['HeII'] += k_HeIIH_cx * nH
            nu['H'] += k_HeIIH_cx * nHeII
            
            # He+ + He -> He + He+ (symmetric charge exchange, energy exchange only)
            k_HeIIHeI_cx = rates.HeIIHeI_charge_exchange(THeII, THeI) * ndamp(nHeII) * CM3_TO_M3
            rate_HeIIHeI_cx = safe_rate(k_HeIIHeI_cx, nHeII, nHeI)
            
            dE['HeII'] += rate_HeIIHeI_cx * (THeI - THeII)
            dE['HeI'] += rate_HeIIHeI_cx * (THeII - THeI)
            # Collision frequencies
            nu['HeII'] += k_HeIIHeI_cx * nHeI
            nu['HeI'] += k_HeIIHeI_cx * nHeII
            
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
            # Collision frequencies
            nu['HeIII'] += k_HeIIIH_cx * nH
            nu['H'] += k_HeIIIH_cx * nHeIII
            
            # He++ + He -> 2 He+ (double charge exchange)
            k_HeIIIHeI_cxa = rates.HeIIIHeI_cx_double(THeIII, THeI) * ndamp(nHeIII) * CM3_TO_M3
            rate_HeIIIHeI_cxa = safe_rate(k_HeIIIHeI_cxa, nHeIII, nHeI)
            
            dn['HeIII'] -= rate_HeIIIHeI_cxa
            dn['HeII'] += 2.0 * rate_HeIIIHeI_cxa  # Creates 2 He+ ions
            dn['HeI'] -= rate_HeIIIHeI_cxa
            
            dE['HeIII'] -= rate_HeIIIHeI_cxa * THeIII
            dE['HeII'] += rate_HeIIIHeI_cxa * (THeIII + THeI)  # Both energies go to He+
            dE['HeI'] -= rate_HeIIIHeI_cxa * THeI
            # Collision frequencies
            nu['HeIII'] += k_HeIIIHeI_cxa * nHeI
            nu['HeI'] += k_HeIIIHeI_cxa * nHeIII
            
            # He++ + He -> He + He++ (symmetric charge exchange, energy exchange only)
            k_HeIIIHeI_cxb = rates.HeIIIHeI_cx_symmetric(THeIII, THeI) * ndamp(nHeIII) * CM3_TO_M3
            rate_HeIIIHeI_cxb = safe_rate(k_HeIIIHeI_cxb, nHeIII, nHeI)
            
            dE['HeIII'] += rate_HeIIIHeI_cxb * (THeI - THeIII)
            dE['HeI'] += rate_HeIIIHeI_cxb * (THeIII - THeI)
            # Collision frequencies
            nu['HeIII'] += k_HeIIIHeI_cxb * nHeI
            nu['HeI'] += k_HeIIIHeI_cxb * nHeIII
    
    # =========================================================================
    # ION-NEUTRAL REACTIONS (bion)
    # Proton impact excitation, ionization, and other ion-neutral collisions
    # =========================================================================
    if include_ion:
        # H+ + H -> H+ + H* excitation channels (energy loss only)
        # REAC311: excitation channel a (10.2 eV)
        k_HiH_exc_a = rates.HiH_excitation_a(THi, TH) * ndamp(nHi_cgs) * CM3_TO_M3
        rate_HiH_exc_a = safe_rate(k_HiH_exc_a, nHi, nH)
        dE['Hi'] -= rate_HiH_exc_a * 10.2
        # Collision frequencies
        nu['Hi'] += k_HiH_exc_a * nH
        nu['H'] += k_HiH_exc_a * nHi
        
        # REAC312: excitation channel b (10.2 eV)
        k_HiH_exc_b = rates.HiH_excitation_b(THi, TH) * ndamp(nHi_cgs) * CM3_TO_M3
        rate_HiH_exc_b = safe_rate(k_HiH_exc_b, nHi, nH)
        dE['Hi'] -= rate_HiH_exc_b * 10.2
        # Collision frequencies
        nu['Hi'] += k_HiH_exc_b * nH
        nu['H'] += k_HiH_exc_b * nHi
        
        # H+ + H -> 2H+ + e (ionization by proton impact)
        # REAC316: 13.6 eV threshold
        k_HiH_ion = rates.HiH_ionization(THi, TH) * ndamp(nH_cgs) * CM3_TO_M3
        rate_HiH_ion = safe_rate(k_HiH_ion, nHi, nH)
        
        dn['H'] -= rate_HiH_ion
        dn['Hi'] += rate_HiH_ion  # Net: creates one new H+
        dn['e'] += rate_HiH_ion
        
        dE['H'] -= rate_HiH_ion * TH
        dE['Hi'] += rate_HiH_ion * TH  # New H+ gets H temperature
        dE['e'] += rate_HiH_ion * 0.1  # Small electron energy from ionization
        # Collision frequencies
        nu['Hi'] += k_HiH_ion * nH
        nu['H'] += k_HiH_ion * nHi
        
        # H2-related ion-neutral reactions
        if include_H2 and nH2 is not None:
            nH2_cgs = nH2 * M3_TO_CM3
            nH2i_cgs = nH2i * M3_TO_CM3 if nH2i is not None else 0.0
            
            # H+ + H2 -> H+ + H2* vibrational excitation (REAC321, 0.1 eV)
            k_HiH2_exc_a = rates.HiH2_excitation_a(THi, TH2) * ndamp(nH2_cgs) * CM3_TO_M3
            rate_HiH2_exc_a = safe_rate(k_HiH2_exc_a, nHi, nH2)
            dE['Hi'] -= rate_HiH2_exc_a * 0.1
            # Collision frequencies
            nu['Hi'] += k_HiH2_exc_a * nH2
            nu['H2'] += k_HiH2_exc_a * nHi
            
            # H+ + H2 -> H+ + H2* electronic excitation (REAC322, 1.0 eV)
            k_HiH2_exc_b = rates.HiH2_excitation_b(THi, TH2) * ndamp(nH2_cgs) * CM3_TO_M3
            rate_HiH2_exc_b = safe_rate(k_HiH2_exc_b, nHi, nH2)
            dE['Hi'] -= rate_HiH2_exc_b * 1.0
            # Collision frequencies
            nu['Hi'] += k_HiH2_exc_b * nH2
            nu['H2'] += k_HiH2_exc_b * nHi
            
            # H+ + H2 -> H+ + H2+ + e ionization (REAC325, 15.4 eV)
            k_HiH2_ion = rates.HiH2_ionization(THi, TH2) * ndamp(nH2_cgs) * CM3_TO_M3
            rate_HiH2_ion = safe_rate(k_HiH2_ion, nHi, nH2)
            
            dn['H2'] -= rate_HiH2_ion
            dn['H2i'] += rate_HiH2_ion
            dn['e'] += rate_HiH2_ion
            
            dE['H2'] -= rate_HiH2_ion * TH2
            dE['H2i'] += rate_HiH2_ion * TH2
            dE['e'] += rate_HiH2_ion * 0.1
            # Collision frequencies
            nu['Hi'] += k_HiH2_ion * nH2
            nu['H2'] += k_HiH2_ion * nHi
            
            if nH2i is not None:
                # H+ + H2+ -> H + H + H+ dissociation (REAC326, 10.5 eV)
                k_HiH2i_diss = rates.HiH2i_dissociation(THi, TH2i) * ndamp(nH2i_cgs) * CM3_TO_M3
                rate_HiH2i_diss = safe_rate(k_HiH2i_diss, nHi, nH2i)
                
                dn['H2i'] -= rate_HiH2i_diss
                dn['H'] += 2.0 * rate_HiH2i_diss  # Creates 2 H atoms
                # H+ is catalyst, no net change
                
                dE['H2i'] -= rate_HiH2i_diss * TH2i
                dE['H'] += rate_HiH2i_diss * TH2i  # Energy to H atoms
                # Collision frequencies
                nu['Hi'] += k_HiH2i_diss * nH2i
                nu['H2i'] += k_HiH2i_diss * nHi
                
                # H2+ + H2 -> H3+ + H (REAC433)
                k_H2iH2_H3i = rates.H2iH2_to_H3i(TH2i, TH2) * ndamp(nH2i_cgs) * CM3_TO_M3
                rate_H2iH2_H3i = safe_rate(k_H2iH2_H3i, nH2i, nH2)
                
                dn['H2i'] -= rate_H2iH2_H3i
                dn['H2'] -= rate_H2iH2_H3i
                dn['H3i'] += rate_H2iH2_H3i
                dn['H'] += rate_H2iH2_H3i
                
                dE['H2i'] -= rate_H2iH2_H3i * TH2i
                dE['H2'] -= rate_H2iH2_H3i * TH2
                dE['H3i'] += rate_H2iH2_H3i * (TH2i + TH2) * 0.75  # 3/4 to H3+
                dE['H'] += rate_H2iH2_H3i * (TH2i + TH2) * 0.25   # 1/4 to H
                # Collision frequencies
                nu['H2i'] += k_H2iH2_H3i * nH2
                nu['H2'] += k_H2iH2_H3i * nH2i
        
        # He-related ion-neutral reactions
        if include_He and nHeI is not None:
            nHeI_cgs = nHeI * M3_TO_CM3
            nHeII_cgs = nHeII * M3_TO_CM3 if nHeII is not None else 0.0
            
            # H+ + He -> H+ + He+ + e ionization (REAC332, 24.58 eV)
            k_HiHeI_ion = rates.HiHeI_ionization(THi, THeI) * ndamp(nHeI_cgs) * CM3_TO_M3
            rate_HiHeI_ion = safe_rate(k_HiHeI_ion, nHi, nHeI)
            
            dn['HeI'] -= rate_HiHeI_ion
            dn['HeII'] += rate_HiHeI_ion
            dn['e'] += rate_HiHeI_ion
            
            dE['HeI'] -= rate_HiHeI_ion * THeI
            dE['HeII'] += rate_HiHeI_ion * THeI
            dE['e'] += rate_HiHeI_ion * 0.1
            # Collision frequencies
            nu['Hi'] += k_HiHeI_ion * nHeI
            nu['HeI'] += k_HiHeI_ion * nHi
            
            # He+ + H2 -> He + H + H+ charge exchange dissociation (REAC523)
            if include_H2 and nH2 is not None and nHeII is not None:
                k_HeIIH2_cxdis = rates.HeIIH2_cx_dissociation(THeII, TH2) * ndamp(nHeII_cgs) * CM3_TO_M3
                rate_HeIIH2_cxdis = safe_rate(k_HeIIH2_cxdis, nHeII, nH2)
                
                dn['HeII'] -= rate_HeIIH2_cxdis
                dn['HeI'] += rate_HeIIH2_cxdis
                dn['H2'] -= rate_HeIIH2_cxdis
                dn['H'] += rate_HeIIH2_cxdis  # One H atom
                dn['Hi'] += rate_HeIIH2_cxdis  # One H+ ion
                # Note: No electron change - charge is conserved (He+ + H2 -> He + H + H+)
                
                dE['HeII'] -= rate_HeIIH2_cxdis * THeII
                dE['HeI'] += rate_HeIIH2_cxdis * THeII
                dE['H2'] -= rate_HeIIH2_cxdis * TH2
                dE['H'] += rate_HeIIH2_cxdis * TH2 * 0.5
                dE['Hi'] += rate_HeIIH2_cxdis * TH2 * 0.5
                # Collision frequencies
                nu['HeII'] += k_HeIIH2_cxdis * nH2
                nu['H2'] += k_HeIIH2_cxdis * nHeII
    
    # =========================================================================
    # ELASTIC COLLISIONS (energy exchange only, no density change)
    # dE = -k * n1 * n2 * L * (T1 - T2)
    # where L = 4*m1*m2/(m1+m2)^2 is the Langevin energy loss factor
    # =========================================================================
    if include_elastic:
        from . import rates as elastic_rates
        
        # H+ + H elastic
        if include_H:
            L_HiH = elastic_rates.Langevin_energy_loss_factor(1.0, 1.0)  # L = 1.0
            k_HiH_el = elastic_rates._elastic_rate_HiH(THi, TH) * CM3_TO_M3
            rate_HiH_el = safe_rate(k_HiH_el, nHi, nH)
            dE['Hi'] -= rate_HiH_el * L_HiH * (THi - TH)
            dE['H'] += rate_HiH_el * L_HiH * (THi - TH)
            # Collision frequencies
            nu['Hi'] += k_HiH_el * nH
            nu['H'] += k_HiH_el * nHi
        
        # H + H self-collision (hard sphere) - contributes to nu for diffusion
        # No energy exchange (same species, same temperature), but important for mean free path
        if include_H:
            k_HH = elastic_rates._elastic_rate_HH(TH) * CM3_TO_M3
            nu['H'] += k_HH * nH
        
        # H2 related elastic collisions
        if include_H2 and nH2 is not None:
            # H + H2 elastic
            L_HH2 = elastic_rates.Langevin_energy_loss_factor(1.0, 2.0)  # L = 0.889
            k_HH2_el = elastic_rates._elastic_rate_HH2(TH, TH2) * CM3_TO_M3
            rate_HH2_el = safe_rate(k_HH2_el, nH, nH2)
            dE['H'] -= rate_HH2_el * L_HH2 * (TH - TH2)
            dE['H2'] += rate_HH2_el * L_HH2 * (TH - TH2)
            # Collision frequencies
            nu['H'] += k_HH2_el * nH2
            nu['H2'] += k_HH2_el * nH
            
            # H+ + H2 elastic
            L_HiH2 = elastic_rates.Langevin_energy_loss_factor(1.0, 2.0)  # L = 0.889
            k_HiH2_el = elastic_rates._elastic_rate_HiH2(THi, TH2) * CM3_TO_M3
            rate_HiH2_el = safe_rate(k_HiH2_el, nHi, nH2)
            dE['Hi'] -= rate_HiH2_el * L_HiH2 * (THi - TH2)
            dE['H2'] += rate_HiH2_el * L_HiH2 * (THi - TH2)
            # Collision frequencies
            nu['Hi'] += k_HiH2_el * nH2
            nu['H2'] += k_HiH2_el * nHi
            
            if nH2i is not None:
                # H2+ + H elastic
                L_H2iH = elastic_rates.Langevin_energy_loss_factor(2.0, 1.0)  # L = 0.889
                k_H2iH_el = elastic_rates._elastic_rate_H2iH(TH2i, TH) * CM3_TO_M3
                rate_H2iH_el = safe_rate(k_H2iH_el, nH2i, nH)
                dE['H2i'] -= rate_H2iH_el * L_H2iH * (TH2i - TH)
                dE['H'] += rate_H2iH_el * L_H2iH * (TH2i - TH)
                # Collision frequencies
                nu['H2i'] += k_H2iH_el * nH
                nu['H'] += k_H2iH_el * nH2i
                
                # H2+ + H2 elastic
                L_H2iH2 = elastic_rates.Langevin_energy_loss_factor(2.0, 2.0)  # L = 1.0
                k_H2iH2_el = elastic_rates._elastic_rate_H2iH2(TH2i, TH2) * CM3_TO_M3
                rate_H2iH2_el = safe_rate(k_H2iH2_el, nH2i, nH2)
                dE['H2i'] -= rate_H2iH2_el * L_H2iH2 * (TH2i - TH2)
                dE['H2'] += rate_H2iH2_el * L_H2iH2 * (TH2i - TH2)
                # Collision frequencies
                nu['H2i'] += k_H2iH2_el * nH2
                nu['H2'] += k_H2iH2_el * nH2i
            
            if nH3i is not None:
                # H3+ + H elastic
                L_H3iH = elastic_rates.Langevin_energy_loss_factor(3.0, 1.0)  # L = 0.75
                k_H3iH_el = elastic_rates._elastic_rate_H3iH(TH3i, TH) * CM3_TO_M3
                rate_H3iH_el = safe_rate(k_H3iH_el, nH3i, nH)
                dE['H3i'] -= rate_H3iH_el * L_H3iH * (TH3i - TH)
                dE['H'] += rate_H3iH_el * L_H3iH * (TH3i - TH)
                # Collision frequencies
                nu['H3i'] += k_H3iH_el * nH
                nu['H'] += k_H3iH_el * nH3i
                
                # H3+ + H2 elastic
                L_H3iH2 = elastic_rates.Langevin_energy_loss_factor(3.0, 2.0)  # L = 0.96
                k_H3iH2_el = elastic_rates._elastic_rate_H3iH2(TH3i, TH2) * CM3_TO_M3
                rate_H3iH2_el = safe_rate(k_H3iH2_el, nH3i, nH2)
                dE['H3i'] -= rate_H3iH2_el * L_H3iH2 * (TH3i - TH2)
                dE['H2'] += rate_H3iH2_el * L_H3iH2 * (TH3i - TH2)
                # Collision frequencies
                nu['H3i'] += k_H3iH2_el * nH2
                nu['H2'] += k_H3iH2_el * nH3i
            
            # H2 + H2 self-collision (hard sphere) - contributes to nu for diffusion
            # No energy exchange (same species, same temperature), but important for mean free path
            k_H2H2 = elastic_rates._elastic_rate_H2H2(TH2) * CM3_TO_M3
            nu['H2'] += k_H2H2 * nH2
        
        # He related elastic collisions
        if include_He and nHeI is not None:
            # H + He elastic
            L_HHe = elastic_rates.Langevin_energy_loss_factor(1.0, 4.0)  # L = 0.64
            k_HHe_el = elastic_rates._elastic_rate_HHe(TH, THeI) * CM3_TO_M3
            rate_HHe_el = safe_rate(k_HHe_el, nH, nHeI)
            dE['H'] -= rate_HHe_el * L_HHe * (TH - THeI)
            dE['HeI'] += rate_HHe_el * L_HHe * (TH - THeI)
            # Collision frequencies
            nu['H'] += k_HHe_el * nHeI
            nu['HeI'] += k_HHe_el * nH
            
            if nHeII is not None:
                # He+ + H elastic
                L_HeIIH = elastic_rates.Langevin_energy_loss_factor(4.0, 1.0)  # L = 0.64
                k_HeIIH_el = elastic_rates._elastic_rate_HeIIH(THeII, TH) * CM3_TO_M3
                rate_HeIIH_el = safe_rate(k_HeIIH_el, nHeII, nH)
                dE['HeII'] -= rate_HeIIH_el * L_HeIIH * (THeII - TH)
                dE['H'] += rate_HeIIH_el * L_HeIIH * (THeII - TH)
                # Collision frequencies
                nu['HeII'] += k_HeIIH_el * nH
                nu['H'] += k_HeIIH_el * nHeII
                
                # H+ + He elastic
                L_HiHe = elastic_rates.Langevin_energy_loss_factor(1.0, 4.0)  # L = 0.64
                k_HiHe_el = elastic_rates._elastic_rate_HiHe(THi, THeI) * CM3_TO_M3
                rate_HiHe_el = safe_rate(k_HiHe_el, nHi, nHeI)
                dE['Hi'] -= rate_HiHe_el * L_HiHe * (THi - THeI)
                dE['HeI'] += rate_HiHe_el * L_HiHe * (THi - THeI)
                # Collision frequencies
                nu['Hi'] += k_HiHe_el * nHeI
                nu['HeI'] += k_HiHe_el * nHi
                
                # He+ + He elastic
                L_HeIIHe = elastic_rates.Langevin_energy_loss_factor(4.0, 4.0)  # L = 1.0
                k_HeIIHe_el = elastic_rates._elastic_rate_HeIIHe(THeII, THeI) * CM3_TO_M3
                rate_HeIIHe_el = safe_rate(k_HeIIHe_el, nHeII, nHeI)
                dE['HeII'] -= rate_HeIIHe_el * L_HeIIHe * (THeII - THeI)
                dE['HeI'] += rate_HeIIHe_el * L_HeIIHe * (THeII - THeI)
                # Collision frequencies
                nu['HeII'] += k_HeIIHe_el * nHeI
                nu['HeI'] += k_HeIIHe_el * nHeII
            
            # HeI + HeI self-collision (hard sphere) - contributes to nu for diffusion
            # No energy exchange (same species, same temperature), but important for mean free path
            k_HeIHeI = elastic_rates._elastic_rate_HeIHeI(THeI) * CM3_TO_M3
            nu['HeI'] += k_HeIHeI * nHeI
            
            # He + H2 elastic (uses HH2 rate as approximation)
            if include_H2 and nH2 is not None:
                L_HeH2 = elastic_rates.Langevin_energy_loss_factor(4.0, 2.0)  # L = 0.889
                k_HeH2_el = elastic_rates._elastic_rate_HH2(THeI, TH2) * CM3_TO_M3  # Approximation
                rate_HeH2_el = safe_rate(k_HeH2_el, nHeI, nH2)
                dE['HeI'] -= rate_HeH2_el * L_HeH2 * (THeI - TH2)
                dE['H2'] += rate_HeH2_el * L_HeH2 * (THeI - TH2)
                # Collision frequencies
                nu['HeI'] += k_HeH2_el * nH2
                nu['H2'] += k_HeH2_el * nHeI
                
                if nHeII is not None:
                    # He+ + H2 elastic
                    L_HeIIH2 = elastic_rates.Langevin_energy_loss_factor(4.0, 2.0)  # L = 0.889
                    k_HeIIH2_el = elastic_rates._elastic_rate_HeIIH2(THeII, TH2) * CM3_TO_M3
                    rate_HeIIH2_el = safe_rate(k_HeIIH2_el, nHeII, nH2)
                    dE['HeII'] -= rate_HeIIH2_el * L_HeIIH2 * (THeII - TH2)
                    dE['H2'] += rate_HeIIH2_el * L_HeIIH2 * (THeII - TH2)
                    # Collision frequencies
                    nu['HeII'] += k_HeIIH2_el * nH2
                    nu['H2'] += k_HeIIH2_el * nHeII
    
    # =========================================================================
    # COULOMB COLLISIONS (energy exchange between charged particles)
    # Critical for thermalization between electrons and ions
    # Q = (T1 - T2) * nu * n1 [eV·m^-3/s]
    # nu_Coulomb also contributes to transport (diffusion)
    # =========================================================================
    if include_coulomb:
        from . import rates as coulomb_rates
        
        # Convert densities to CGS for Coulomb collision calculations
        ne_cgs = ne * M3_TO_CM3
        nHi_cgs = nHi * M3_TO_CM3
        
        # Only compute Coulomb collisions if electron density is significant
        # (ne >= 10 cm^-3 threshold from C++ code)
        mask = ne_cgs >= 10.0
        
        if np.any(mask):
            # Electron-H+ Coulomb collisions
            nu_eHi = np.zeros_like(ne)
            nu_eHi[mask] = coulomb_rates.nu_ei(ne_cgs[mask], Te[mask], 
                                               nHi_cgs[mask], THi[mask], 
                                               Z=1.0, mu=1.0)
            QeHi = (Te - THi) * nu_eHi * ne
            # Convert from eV·cm^-3/s to eV·m^-3/s
            QeHi *= CM3_TO_M3
            dE['e'] -= QeHi
            dE['Hi'] += QeHi
            # Collision frequencies (Coulomb)
            nu['e'] += nu_eHi
            nu['Hi'] += nu_eHi
            
            # H2 molecular ions
            if nH2i is not None:
                nH2i_cgs = nH2i * M3_TO_CM3
                
                # Electron-H2+ Coulomb collisions
                nu_eH2i = np.zeros_like(ne)
                nu_eH2i[mask] = coulomb_rates.nu_ei(ne_cgs[mask], Te[mask],
                                                    nH2i_cgs[mask], TH2i[mask],
                                                    Z=1.0, mu=2.0)
                QeH2i = (Te - TH2i) * nu_eH2i * ne
                QeH2i *= CM3_TO_M3
                dE['e'] -= QeH2i
                dE['H2i'] += QeH2i
                nu['e'] += nu_eH2i
                nu['H2i'] += nu_eH2i
                
                # H+ - H2+ Coulomb collisions
                nu_HiH2i = np.zeros_like(ne)
                nu_HiH2i[mask] = coulomb_rates.nu_ii(nHi_cgs[mask], THi[mask], 1.0, 1.0,
                                                     nH2i_cgs[mask], TH2i[mask], 1.0, 2.0)
                Q_HiH2i = (THi - TH2i) * nu_HiH2i * nHi
                Q_HiH2i *= CM3_TO_M3
                dE['Hi'] -= Q_HiH2i
                dE['H2i'] += Q_HiH2i
                nu['Hi'] += nu_HiH2i
                nu['H2i'] += nu_HiH2i
            
            if nH3i is not None:
                nH3i_cgs = nH3i * M3_TO_CM3
                
                # Electron-H3+ Coulomb collisions
                nu_eH3i = np.zeros_like(ne)
                nu_eH3i[mask] = coulomb_rates.nu_ei(ne_cgs[mask], Te[mask],
                                                    nH3i_cgs[mask], TH3i[mask],
                                                    Z=1.0, mu=3.0)
                QeH3i = (Te - TH3i) * nu_eH3i * ne
                QeH3i *= CM3_TO_M3
                dE['e'] -= QeH3i
                dE['H3i'] += QeH3i
                nu['e'] += nu_eH3i
                nu['H3i'] += nu_eH3i
                
                # H+ - H3+ Coulomb collisions
                nu_HiH3i = np.zeros_like(ne)
                nu_HiH3i[mask] = coulomb_rates.nu_ii(nHi_cgs[mask], THi[mask], 1.0, 1.0,
                                                     nH3i_cgs[mask], TH3i[mask], 1.0, 3.0)
                Q_HiH3i = (THi - TH3i) * nu_HiH3i * nHi
                Q_HiH3i *= CM3_TO_M3
                dE['Hi'] -= Q_HiH3i
                dE['H3i'] += Q_HiH3i
                nu['Hi'] += nu_HiH3i
                nu['H3i'] += nu_HiH3i
                
                if nH2i is not None:
                    # H2+ - H3+ Coulomb collisions
                    nu_H2iH3i = np.zeros_like(ne)
                    nu_H2iH3i[mask] = coulomb_rates.nu_ii(nH2i_cgs[mask], TH2i[mask], 1.0, 2.0,
                                                          nH3i_cgs[mask], TH3i[mask], 1.0, 3.0)
                    Q_H2iH3i = (TH2i - TH3i) * nu_H2iH3i * nH2i
                    Q_H2iH3i *= CM3_TO_M3
                    dE['H2i'] -= Q_H2iH3i
                    dE['H3i'] += Q_H2iH3i
                    nu['H2i'] += nu_H2iH3i
                    nu['H3i'] += nu_H2iH3i
            
            # Helium ions
            if nHeII is not None:
                nHeII_cgs = nHeII * M3_TO_CM3
                
                # Electron-He+ Coulomb collisions (Z=1)
                nu_eHeII = np.zeros_like(ne)
                nu_eHeII[mask] = coulomb_rates.nu_ei(ne_cgs[mask], Te[mask],
                                                     nHeII_cgs[mask], THeII[mask],
                                                     Z=1.0, mu=4.0)
                QeHeII = (Te - THeII) * nu_eHeII * ne
                QeHeII *= CM3_TO_M3
                dE['e'] -= QeHeII
                dE['HeII'] += QeHeII
                nu['e'] += nu_eHeII
                nu['HeII'] += nu_eHeII
                
                # H+ - He+ Coulomb collisions
                nu_HiHeII = np.zeros_like(ne)
                nu_HiHeII[mask] = coulomb_rates.nu_ii(nHi_cgs[mask], THi[mask], 1.0, 1.0,
                                                      nHeII_cgs[mask], THeII[mask], 1.0, 4.0)
                Q_HiHeII = (THi - THeII) * nu_HiHeII * nHi
                Q_HiHeII *= CM3_TO_M3
                dE['Hi'] -= Q_HiHeII
                dE['HeII'] += Q_HiHeII
                nu['Hi'] += nu_HiHeII
                nu['HeII'] += nu_HiHeII
                
                if nH2i is not None:
                    # H2+ - He+ Coulomb collisions
                    nu_H2iHeII = np.zeros_like(ne)
                    nu_H2iHeII[mask] = coulomb_rates.nu_ii(nH2i_cgs[mask], TH2i[mask], 1.0, 2.0,
                                                           nHeII_cgs[mask], THeII[mask], 1.0, 4.0)
                    Q_H2iHeII = (TH2i - THeII) * nu_H2iHeII * nH2i
                    Q_H2iHeII *= CM3_TO_M3
                    dE['H2i'] -= Q_H2iHeII
                    dE['HeII'] += Q_H2iHeII
                    nu['H2i'] += nu_H2iHeII
                    nu['HeII'] += nu_H2iHeII
                
                if nH3i is not None:
                    # H3+ - He+ Coulomb collisions
                    nu_H3iHeII = np.zeros_like(ne)
                    nu_H3iHeII[mask] = coulomb_rates.nu_ii(nH3i_cgs[mask], TH3i[mask], 1.0, 3.0,
                                                           nHeII_cgs[mask], THeII[mask], 1.0, 4.0)
                    Q_H3iHeII = (TH3i - THeII) * nu_H3iHeII * nH3i
                    Q_H3iHeII *= CM3_TO_M3
                    dE['H3i'] -= Q_H3iHeII
                    dE['HeII'] += Q_H3iHeII
                    nu['H3i'] += nu_H3iHeII
                    nu['HeII'] += nu_H3iHeII
            
            if nHeIII is not None:
                nHeIII_cgs = nHeIII * M3_TO_CM3
                
                # Electron-He++ Coulomb collisions (Z=2)
                nu_eHeIII = np.zeros_like(ne)
                nu_eHeIII[mask] = coulomb_rates.nu_ei(ne_cgs[mask], Te[mask],
                                                      nHeIII_cgs[mask], THeIII[mask],
                                                      Z=2.0, mu=4.0)
                QeHeIII = (Te - THeIII) * nu_eHeIII * ne
                QeHeIII *= CM3_TO_M3
                dE['e'] -= QeHeIII
                dE['HeIII'] += QeHeIII
                nu['e'] += nu_eHeIII
                nu['HeIII'] += nu_eHeIII
                
                # H+ - He++ Coulomb collisions
                nu_HiHeIII = np.zeros_like(ne)
                nu_HiHeIII[mask] = coulomb_rates.nu_ii(nHi_cgs[mask], THi[mask], 1.0, 1.0,
                                                       nHeIII_cgs[mask], THeIII[mask], 2.0, 4.0)
                Q_HiHeIII = (THi - THeIII) * nu_HiHeIII * nHi
                Q_HiHeIII *= CM3_TO_M3
                dE['Hi'] -= Q_HiHeIII
                dE['HeIII'] += Q_HiHeIII
                nu['Hi'] += nu_HiHeIII
                nu['HeIII'] += nu_HiHeIII
                
                if nH2i is not None:
                    # H2+ - He++ Coulomb collisions
                    nu_H2iHeIII = np.zeros_like(ne)
                    nu_H2iHeIII[mask] = coulomb_rates.nu_ii(nH2i_cgs[mask], TH2i[mask], 1.0, 2.0,
                                                            nHeIII_cgs[mask], THeIII[mask], 2.0, 4.0)
                    Q_H2iHeIII = (TH2i - THeIII) * nu_H2iHeIII * nH2i
                    Q_H2iHeIII *= CM3_TO_M3
                    dE['H2i'] -= Q_H2iHeIII
                    dE['HeIII'] += Q_H2iHeIII
                    nu['H2i'] += nu_H2iHeIII
                    nu['HeIII'] += nu_H2iHeIII
                
                if nH3i is not None:
                    # H3+ - He++ Coulomb collisions
                    nu_H3iHeIII = np.zeros_like(ne)
                    nu_H3iHeIII[mask] = coulomb_rates.nu_ii(nH3i_cgs[mask], TH3i[mask], 1.0, 3.0,
                                                            nHeIII_cgs[mask], THeIII[mask], 2.0, 4.0)
                    Q_H3iHeIII = (TH3i - THeIII) * nu_H3iHeIII * nH3i
                    Q_H3iHeIII *= CM3_TO_M3
                    dE['H3i'] -= Q_H3iHeIII
                    dE['HeIII'] += Q_H3iHeIII
                    nu['H3i'] += nu_H3iHeIII
                    nu['HeIII'] += nu_H3iHeIII
                
                if nHeII is not None:
                    # He+ - He++ Coulomb collisions
                    nu_HeIIHeIII = np.zeros_like(ne)
                    nu_HeIIHeIII[mask] = coulomb_rates.nu_ii(nHeII_cgs[mask], THeII[mask], 1.0, 4.0,
                                                             nHeIII_cgs[mask], THeIII[mask], 2.0, 4.0)
                    Q_HeIIHeIII = (THeII - THeIII) * nu_HeIIHeIII * nHeII
                    Q_HeIIHeIII *= CM3_TO_M3
                    dE['HeII'] -= Q_HeIIHeIII
                    dE['HeIII'] += Q_HeIIHeIII
                    nu['HeII'] += nu_HeIIHeIII
                    nu['HeIII'] += nu_HeIIHeIII
    
    return dn, dE, nu


def compute_sources_from_state(state, params: dict) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], Dict[str, np.ndarray]]:
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
    nu : Dict[str, np.ndarray]
        Collision frequencies [1/s].
    """
    # Extract density and temperature arrays
    ne = state.electrons.n.x.array.copy()
    Te = state.electrons.T.copy()
    
    nH = state.species['H'].n.x.array.copy() if 'H' in state.species else np.zeros_like(ne)
    TH = state.species['H'].T.copy() if 'H' in state.species else np.ones_like(ne)
    
    nHi = state.species['Hi'].n.x.array.copy() if 'Hi' in state.species else np.zeros_like(ne)
    THi = state.species['Hi'].T.copy() if 'Hi' in state.species else Te.copy()
    
    # Molecular hydrogen species (optional)
    nH2 = state.species['H2'].n.x.array.copy() if 'H2' in state.species else None
    TH2 = state.species['H2'].T.copy() if 'H2' in state.species else None
    
    nH2i = state.species['H2i'].n.x.array.copy() if 'H2i' in state.species else None
    TH2i = state.species['H2i'].T.copy() if 'H2i' in state.species else None
    
    nH3i = state.species['H3i'].n.x.array.copy() if 'H3i' in state.species else None
    TH3i = state.species['H3i'].T.copy() if 'H3i' in state.species else None
    
    # Helium species (optional)
    nHeI = state.species['HeI'].n.x.array.copy() if 'HeI' in state.species else None
    THeI = state.species['HeI'].T.copy() if 'HeI' in state.species else None
    
    nHeII = state.species['HeII'].n.x.array.copy() if 'HeII' in state.species else None
    THeII = state.species['HeII'].T.copy() if 'HeII' in state.species else None
    
    nHeIII = state.species['HeIII'].n.x.array.copy() if 'HeIII' in state.species else None
    THeIII = state.species['HeIII'].T.copy() if 'HeIII' in state.species else None
    
    # Get physics flags from params
    include_H = params.get('bH', True)
    include_H2 = params.get('bH2', 'H2' in state.species)
    include_He = params.get('bHe', 'HeI' in state.species)
    include_ion = params.get('bion', True)
    include_cx = params.get('bcx', True)
    include_elastic = params.get('belas', True)
    include_coulomb = params.get('bcoulomb', True)
    use_ADAS = params.get('bADAS', True)
    
    return compute_collision_sources(
        ne, Te, nH, TH, nHi, THi,
        nH2, TH2, nH2i, TH2i, nH3i, TH3i,
        nHeI, THeI, nHeII, THeII, nHeIII, THeIII,
        use_ADAS=use_ADAS,
        include_H=include_H,
        include_H2=include_H2,
        include_He=include_He,
        include_ion=include_ion,
        include_cx=include_cx,
        include_elastic=include_elastic,
        include_coulomb=include_coulomb
    )
