"""
Reaction rate coefficients for plasma collisions.

Implements rate coefficients from HYDHEL/IAEA databases for:
- Electron-impact ionization, excitation, recombination
- Charge exchange reactions
- Ion-impact reactions

Rate coefficients are in cm³/s (CGS units, matching C++ Tomator1D).

For He reactions, uses exact ADAS 2D interpolation tables from C++ reactionrates.cpp.
"""

import numpy as np
from typing import Tuple

# Import ADAS tables for He reactions
from .rates_adas import (
    HeI_ionization_ADAS, HeII_ionization_ADAS,
    HeII_recombination_ADAS, HeIII_recombination_ADAS,
    He_cooling_rate_ADAS
)

# Physical constants
QE = 1.60218e-19   # Elementary charge [C]
ME = 9.10938e-31   # Electron mass [kg]
MI = 1.67262e-27   # Proton mass [kg]


def _safe_exp(x: np.ndarray, max_val: float = 500.0) -> np.ndarray:
    """Safe exponential to avoid overflow."""
    return np.exp(np.clip(x, -max_val, max_val))


def _safe_log(x: np.ndarray, min_val: float = 1e-30) -> np.ndarray:
    """Safe logarithm to avoid -inf."""
    return np.log(np.maximum(x, min_val))


class ReactionRates:
    """
    Reaction rate coefficient calculations from HYDHEL database.
    
    All rates are returned in cm³/s (CGS units).
    Temperature inputs are in eV.
    """
    
    # -------------------------------------------------------------------
    # Electron-Hydrogen Atom Reactions (bH)
    # -------------------------------------------------------------------
    
    @staticmethod
    def H_excitation(Te: np.ndarray) -> np.ndarray:
        """
        Electron-impact excitation of atomic hydrogen.
        e + H -> e + H*  (Reaction 2.1.1-2.1.4)
        
        Parameters
        ----------
        Te : np.ndarray
            Electron temperature [eV].
            
        Returns
        -------
        k : np.ndarray
            Rate coefficient [cm³/s].
        """
        k = np.zeros_like(Te)
        mask = Te > 0.6
        Te_safe = np.maximum(Te[mask], 0.6)
        
        # From collisions.cpp line 185
        k[mask] = (9.70346e-8 * np.power(10.2 / Te_safe, 0.92457) * 
                   _safe_exp(-10.2 / Te_safe) / (0.01351 + 10.2 / Te_safe))
        return k
    
    @staticmethod
    def H_ionization(Te: np.ndarray) -> np.ndarray:
        """
        Electron-impact ionization of atomic hydrogen.
        e + H -> e + H+ + e  (Reaction 2.1.5-2.1.7)
        
        Parameters
        ----------
        Te : np.ndarray
            Electron temperature [eV].
            
        Returns
        -------
        k : np.ndarray
            Rate coefficient [cm³/s].
        """
        k = np.zeros_like(Te)
        mask = Te > 0.6
        Te_safe = np.maximum(Te[mask], 0.6)
        
        # From collisions.cpp line 193
        k[mask] = (2.91e-8 * np.power(13.6 / Te_safe, 0.39) * 
                   _safe_exp(-13.6 / Te_safe) / (0.232 + 13.6 / Te_safe))
        return k
    
    @staticmethod
    def H_3body_recombination(Te: np.ndarray, k_ion: np.ndarray) -> np.ndarray:
        """
        Three-body recombination of H+.
        e + e + H+ -> e + H  (reverse of ionization)
        
        Parameters
        ----------
        Te : np.ndarray
            Electron temperature [eV].
        k_ion : np.ndarray
            Ionization rate coefficient [cm³/s].
            
        Returns
        -------
        k : np.ndarray
            Rate coefficient [cm⁶/s].
        """
        # From collisions.cpp line 205: detailed balance
        # 170.8160 = 4 * pi * 13.6
        k = 1.4804e-25 * np.power(170.8160 / np.maximum(Te, 0.1), 1.5) * _safe_exp(13.6 / np.maximum(Te, 0.1)) * k_ion
        return k
    
    @staticmethod
    def H_radiative_recombination(Te: np.ndarray) -> np.ndarray:
        """
        Radiative recombination of H+.
        e + H+ -> H + photon  (Reaction 2.1.8a)
        
        Parameters
        ----------
        Te : np.ndarray
            Electron temperature [eV].
            
        Returns
        -------
        k : np.ndarray
            Rate coefficient [cm³/s].
        """
        k = np.zeros_like(Te)
        mask = Te < 1000.0
        Te_safe = np.maximum(Te[mask], 0.01)
        
        # From collisions.cpp line 220
        sqrt_Te_ratio = np.sqrt(Te_safe / 2.713e-4)
        k[mask] = (7.982e-11 / 
                   (sqrt_Te_ratio * 
                    np.power(1 + sqrt_Te_ratio, 1.0 - 0.7480) * 
                    np.power(1 + np.sqrt(Te_safe / 60.631), 1.0 + 0.7480)))
        return k
    
    # -------------------------------------------------------------------
    # Helium Reactions (bHe)
    # -------------------------------------------------------------------
    
    @staticmethod
    def HeI_ionization(Te: np.ndarray, ne_cgs: np.ndarray = None) -> np.ndarray:
        """
        Electron-impact ionization of neutral helium.
        e + He -> e + He+ + e  (Reaction 2.3.9-2.3.12)
        
        Uses ADAS 2D tables from C++ Tomator1D (Te, ne dependent).
        
        Parameters
        ----------
        Te : np.ndarray
            Electron temperature [eV].
        ne_cgs : np.ndarray, optional
            Electron density [cm^-3]. If None, uses 1e11 cm^-3.
            
        Returns
        -------
        k : np.ndarray
            Rate coefficient [cm³/s].
        """
        if ne_cgs is None:
            ne_cgs = np.full_like(Te, 1e11)
        return HeI_ionization_ADAS(Te, ne_cgs)
    
    @staticmethod
    def HeII_ionization(Te: np.ndarray, ne_cgs: np.ndarray = None) -> np.ndarray:
        """
        Electron-impact ionization of He+.
        e + He+ -> e + He++ + e  (Reaction 2.3.19)
        
        Uses ADAS 2D tables from C++ Tomator1D (Te, ne dependent).
        
        Parameters
        ----------
        Te : np.ndarray
            Electron temperature [eV].
        ne_cgs : np.ndarray, optional
            Electron density [cm^-3]. If None, uses 1e11 cm^-3.
            
        Returns
        -------
        k : np.ndarray
            Rate coefficient [cm³/s].
        """
        if ne_cgs is None:
            ne_cgs = np.full_like(Te, 1e11)
        return HeII_ionization_ADAS(Te, ne_cgs)
    
    @staticmethod
    def HeII_recombination(Te: np.ndarray, ne_cgs: np.ndarray = None) -> np.ndarray:
        """
        Radiative recombination of He+.
        e + He+ -> He + photon  (Reaction 2.3.13)
        
        Uses ADAS 2D tables from C++ Tomator1D (Te, ne dependent).
        
        Parameters
        ----------
        Te : np.ndarray
            Electron temperature [eV].
        ne_cgs : np.ndarray, optional
            Electron density [cm^-3]. If None, uses 1e11 cm^-3.
            
        Returns
        -------
        k : np.ndarray
            Rate coefficient [cm³/s].
        """
        if ne_cgs is None:
            ne_cgs = np.full_like(Te, 1e11)
        return HeII_recombination_ADAS(Te, ne_cgs)
    
    @staticmethod
    def HeIII_recombination(Te: np.ndarray, ne_cgs: np.ndarray = None) -> np.ndarray:
        """
        Radiative recombination of He++.
        e + He++ -> He+ + photon
        
        Uses ADAS 2D tables from C++ Tomator1D (Te, ne dependent).
        
        Parameters
        ----------
        Te : np.ndarray
            Electron temperature [eV].
        ne_cgs : np.ndarray, optional
            Electron density [cm^-3]. If None, uses 1e11 cm^-3.
            
        Returns
        -------
        k : np.ndarray
            Rate coefficient [cm³/s].
        """
        if ne_cgs is None:
            ne_cgs = np.full_like(Te, 1e11)
        return HeIII_recombination_ADAS(Te, ne_cgs)
    
    @staticmethod
    def He_cooling_rate(Te: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Helium cooling rates from ADAS/SOLPS-ITER.
        
        Parameters
        ----------
        Te : np.ndarray
            Electron temperature [eV].
            
        Returns
        -------
        L_HeI : np.ndarray
            Cooling rate from HeI processes [eV*cm³/s].
        L_HeII : np.ndarray
            Cooling rate from HeII processes [eV*cm³/s].
        L_HeIII : np.ndarray
            Cooling rate from HeIII processes [eV*cm³/s].
        """
        # From collisions.cpp lines 479-484 (ADAS cooling tables)
        MTe = np.array([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 1.9, 2.0, 
                        3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0,
                        25.0, 30.0, 40.0, 50.0, 70.0, 100.0, 200.0, 400.0, 
                        600.0, 800.0, 1000.0])
        
        MkHeI = np.array([1.985e-21, 7.830e-21, 2.874e-20, 9.070e-20, 2.055e-19,
                          6.043e-19, 1.176e-18, 4.425e-18, 8.611e-18, 1.299e-17,
                          3.479e-16, 1.672e-15, 5.449e-15, 1.074e-14, 1.836e-14,
                          2.656e-14, 3.659e-14, 4.847e-14, 1.115e-13, 1.784e-13,
                          2.360e-13, 2.931e-13, 3.837e-13, 4.687e-13, 5.793e-13,
                          6.825e-13, 8.126e-13, 8.569e-13, 8.748e-13, 8.847e-13,
                          8.915e-13])  # [eV*m³/s]
        
        MkHeII = np.array([1.058e-19, 1.109e-19, 1.158e-19, 1.206e-19, 1.252e-19,
                           1.296e-19, 1.340e-19, 1.428e-19, 1.478e-19, 1.530e-19,
                           1.171e-18, 1.735e-17, 1.526e-16, 5.449e-16, 1.423e-15,
                           2.707e-15, 4.742e-15, 7.710e-15, 2.857e-14, 5.694e-14,
                           8.292e-14, 1.110e-13, 1.546e-13, 1.982e-13, 2.529e-13,
                           3.058e-13, 3.793e-13, 4.065e-13, 4.072e-13, 3.991e-13,
                           3.897e-13])  # [eV*m³/s]
        
        # Interpolate - keep in eV*m³/s for SI calculations
        L_HeI = np.interp(Te, MTe, MkHeI)  # [eV*m³/s]
        L_HeII = np.interp(Te, MTe, MkHeII)  # [eV*m³/s]
        L_HeIII = np.zeros_like(Te)  # HeIII cooling is negligible
        
        return L_HeI, L_HeII, L_HeIII
    
    @staticmethod
    def He_cooling_rate_IAEA(Te: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Helium cooling rates from IAEA formulas.
        
        Parameters
        ----------
        Te : np.ndarray
            Electron temperature [eV].
            
        Returns
        -------
        L_HeI : np.ndarray
            Cooling rate from HeI processes [W*cm³] (needs /qe for eV*cm³/s).
        L_HeII : np.ndarray
            Cooling rate from HeII processes [W*cm³].
        L_HeIII : np.ndarray
            Cooling rate from HeIII processes [W*cm³].
        """
        # From collisions.cpp lines 500-520 (IAEA cooling)
        Te_safe = np.maximum(Te, 0.1)
        Tx = Te_safe / 1000.0  # Temperature in keV
        
        # HeI coefficients
        HeIa1, HeIa2, HeIa3 = 0.6623e3, 0.9476e-1, 0.7456
        HeIa4, HeIa5, HeIa6 = -0.2592, 3.8098, 0.4026
        
        # HeII coefficients
        HeIIa1, HeIIa2, HeIIa3 = 0.3476e3, 0.1214, 0.7974
        HeIIa4, HeIIa5, HeIIa6 = 0.4819, 1.4066, -0.3639e-2
        
        # Cooling rates [W*cm³]
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            L_HeI = (HeIa1 * _safe_exp(-HeIa2 / np.power(Tx, HeIa3)) / 
                     (np.power(Tx, HeIa4) + HeIa5 * np.power(Tx, HeIa6))) * 1e-33 * 1e6
            
            L_HeII = (HeIIa1 * _safe_exp(-HeIIa2 / np.power(Tx, HeIIa3)) / 
                      (np.power(Tx, HeIIa4) + HeIIa5 * np.power(Tx, HeIIa6))) * 1e-33 * 1e6
        
        L_HeIII = np.zeros_like(Te)
        
        # Handle edge cases
        L_HeI = np.nan_to_num(L_HeI, nan=0.0, posinf=0.0, neginf=0.0)
        L_HeII = np.nan_to_num(L_HeII, nan=0.0, posinf=0.0, neginf=0.0)
        
        return L_HeI, L_HeII, L_HeIII
    
    # -------------------------------------------------------------------
    # Charge Exchange Reactions (bcx)
    # -------------------------------------------------------------------
    
    @staticmethod
    def HiH_charge_exchange(THi: np.ndarray, TH: np.ndarray) -> np.ndarray:
        """
        Charge exchange between H+ and H.
        H+ + H -> H + H+  (Reaction 3.1.8-3.1.11)
        
        Only affects energy transfer, not particle numbers.
        
        Parameters
        ----------
        THi : np.ndarray
            H+ temperature [eV].
        TH : np.ndarray
            H temperature [eV].
            
        Returns
        -------
        k : np.ndarray
            Rate coefficient [cm³/s].
        """
        T_avg = (THi + TH) / 2.0
        k = np.zeros_like(T_avg)
        mask = T_avg > 0.1
        
        # From collisions.cpp line 577
        k[mask] = 7.829e-9 * np.power(T_avg[mask], 0.41)
        return k
    
    @staticmethod
    def HeIIH_charge_exchange(THeII: np.ndarray, TH: np.ndarray) -> np.ndarray:
        """
        Charge exchange: He+ + H -> He + H+
        
        Parameters
        ----------
        THeII : np.ndarray
            He+ temperature [eV].
        TH : np.ndarray
            H temperature [eV].
            
        Returns
        -------
        k : np.ndarray
            Rate coefficient [cm³/s].
        """
        T_avg = (THeII + TH) / 2.0
        T_safe = np.maximum(T_avg, 0.1)
        
        # Simplified fit - actual uses RRCX lookup table
        k = 1.0e-9 * np.power(T_safe / 10.0, 0.3)
        return k
    
    @staticmethod
    def HeIIIH_charge_exchange(THeIII: np.ndarray, TH: np.ndarray) -> np.ndarray:
        """
        Charge exchange: He++ + H -> He+ + H+
        
        Parameters
        ----------
        THeIII : np.ndarray
            He++ temperature [eV].
        TH : np.ndarray
            H temperature [eV].
            
        Returns
        -------
        k : np.ndarray
            Rate coefficient [cm³/s].
        """
        T_avg = (THeIII + TH) / 2.0
        T_safe = np.maximum(T_avg, 0.1)
        
        # Simplified fit
        k = 2.0e-9 * np.power(T_safe / 10.0, 0.3)
        return k
    
    # -------------------------------------------------------------------
    # Elastic Collisions
    # -------------------------------------------------------------------
    
    @staticmethod
    def electron_H2_elastic(ne: np.ndarray, Te: np.ndarray) -> np.ndarray:
        """
        Electron-H2 elastic collision rate.
        
        Parameters
        ----------
        ne : np.ndarray
            Electron density [cm^-3].
        Te : np.ndarray
            Electron temperature [eV].
            
        Returns
        -------
        k : np.ndarray
            Rate coefficient [cm³/s].
        """
        Te_safe = np.maximum(Te, 0.1)
        # Simplified momentum transfer rate
        k = 1.0e-8 * np.power(Te_safe, 0.5)
        return k
    
    @staticmethod
    def Langevin_energy_loss_factor(m1: float, m2: float) -> float:
        """
        Langevin energy loss parameter for elastic collisions.
        
        L = 4 * m1 * m2 / (m1 + m2)²
        
        Parameters
        ----------
        m1 : float
            Mass of species 1 [amu].
        m2 : float
            Mass of species 2 [amu].
            
        Returns
        -------
        L : float
            Energy loss parameter (dimensionless).
        """
        return 4.0 * m1 * m2 / (m1 + m2)**2


def ndamp(n_cgs: np.ndarray, nevac: float = 1.0) -> np.ndarray:
    """
    Density damping factor from C++ Tomator1D (functions.cpp line 392).
    
    ndamp(n) = (1 + nevac/n)^(-0.66)
    
    This limits reactions at very low densities (approaching vacuum).
    
    Parameters
    ----------
    n_cgs : np.ndarray
        Density [cm^-3].
    nevac : float
        Vacuum density reference [cm^-3]. Default is 1.0 from C++ simparam.cpp.
        
    Returns
    -------
    factor : np.ndarray
        Damping factor (0 to 1).
    """
    n_safe = np.maximum(n_cgs, 1e-10)  # Prevent division by zero
    return np.power(1.0 + nevac / n_safe, -0.66)
