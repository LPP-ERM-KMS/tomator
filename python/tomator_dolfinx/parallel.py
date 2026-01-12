"""
Parallel transport losses to limiters in the SOL region.

In the scrape-off layer (R < lHFS or R > lLFS), ions are lost along field 
lines to limiters with a characteristic parallel loss time. Lost ions 
recycle as neutrals at wall temperature.

This module implements the limiters() function from C++ Tomator1D.

Key functions:
- compute_parallel_loss_rates(): Returns k_n, k_E arrays (1/τ) for implicit LHS treatment
- compute_limiter_losses(): Legacy function returning explicit dn/dE source terms
- compute_bpol_losses(): Legacy function returning explicit dn/dE source terms
"""

from typing import Dict, Tuple
import numpy as np


# Unit conversion: m^-3 to cm^-3
M3_TO_CM3 = 1e-6


def compute_parallel_loss_rates(
    R_positions: np.ndarray,
    lHFS: float,
    lLFS: float,
    D_perp: np.ndarray,
    lambda_n: float,
    lambda_E: float,
    Br: np.ndarray = None,
    Bv: float = None,
    b: float = None,
    ne: np.ndarray = None,
    Te: np.ndarray = None,
    nHi: np.ndarray = None,
    THi: np.ndarray = None,
    nuHi: np.ndarray = None,
    nH2i: np.ndarray = None,
    TH2i: np.ndarray = None,
    nuH2i: np.ndarray = None,
    nH3i: np.ndarray = None,
    TH3i: np.ndarray = None,
    nuH3i: np.ndarray = None,
    nHeII: np.ndarray = None,
    THeII: np.ndarray = None,
    nuHeII: np.ndarray = None,
    nHeIII: np.ndarray = None,
    THeIII: np.ndarray = None,
    nuHeIII: np.ndarray = None,
    Dfsave: float = 1.0,
    D_ion: np.ndarray = None,
    diffusion_model: str = 'fixed',
    gEd: float = 1.0
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute parallel loss rates k_n = 1/τ_n and k_E = 1/τ_E for implicit treatment.
    
    Combines limiter losses (SOL region) and bpol losses (vertical diffusion).
    These rates are the SAME for all charged species, allowing implicit treatment
    on the LHS of the transport equation: dn/dt + k_n * n = ...
    
    The implicit formulation is unconditionally stable for any timestep.
    
    Parameters
    ----------
    R_positions : np.ndarray
        Radial positions [m].
    lHFS : float
        HFS limiter position [m].
    lLFS : float
        LFS limiter position [m].
    D_perp : np.ndarray
        Perpendicular diffusion coefficient [m²/s].
    lambda_n : float
        Density decay length for limiter losses [m].
    lambda_E : float
        Energy decay length for limiter losses [m].
    Br : np.ndarray, optional
        Local toroidal magnetic field [T] for bpol losses.
    Bv : float, optional
        Vertical magnetic field [T] for bpol losses.
    b : float, optional
        Vertical plasma extent [cm] for bpol losses.
    ne, Te : np.ndarray, optional
        Electron density [m^-3] and temperature [eV] for gyrogeom model.
    nHi, THi, nuHi : np.ndarray, optional
        H+ quantities for gyrogeom model.
    ... (similar for other species)
    Dfsave : float
        Diffusion scaling factor for gyrogeom model. Default 1.0.
    D_ion : np.ndarray, optional
        Pre-computed ion diffusion coefficient [m²/s] for fixed/bohm models.
    diffusion_model : str
        'fixed', 'bohm', or 'gyrogeom'. Default 'fixed'.
    gEd : float
        Energy loss rate factor. Default 1.0.
        
    Returns
    -------
    k_n : np.ndarray
        Density loss rate 1/τ_n [s^-1]. Same for all charged species.
    k_E : np.ndarray
        Energy loss rate 1/τ_E [s^-1]. Same for all charged species.
    """
    n_points = len(R_positions)
    
    # Initialize loss rates (will sum contributions from limiter + bpol)
    k_n = np.zeros(n_points)
    k_E = np.zeros(n_points)
    
    # ================================================================
    # LIMITER LOSSES (SOL region only)
    # τ_n = λ_n² / D_⊥,  τ_E = λ_E² / D_⊥
    # ================================================================
    
    # Find SOL mask: points outside the confined region
    # Exclude boundary points (handled by decay length BCs)
    sol_mask = (R_positions < lHFS) | (R_positions > lLFS)
    sol_mask[0] = False   # HFS boundary
    sol_mask[-1] = False  # LFS boundary
    
    if np.any(sol_mask):
        D_safe = np.maximum(D_perp, 1e-10)
        tau_n_lim = lambda_n**2 / D_safe
        tau_E_lim = lambda_E**2 / D_safe
        
        tau_n_lim = np.maximum(tau_n_lim, 1e-10)
        tau_E_lim = np.maximum(tau_E_lim, 1e-10)
        
        k_n[sol_mask] += 1.0 / tau_n_lim[sol_mask]
        k_E[sol_mask] += 1.0 / tau_E_lim[sol_mask]
    
    # ================================================================
    # BPOL LOSSES (vertical diffusion, entire domain)
    # τ = b² / (2 * Dv)
    # ================================================================
    
    if Br is not None and Bv is not None and b is not None:
        # Compute vertical diffusion coefficient
        Br_safe = np.maximum(np.abs(Br), 1e-6)
        
        if diffusion_model.lower() in ('fixed', 'bohm') and D_ion is not None:
            # Use pre-computed diffusion coefficient
            Dv = D_ion * 1e4  # [cm²/s]
        elif diffusion_model.lower() == 'gyrogeom' and ne is not None:
            # Compute gyrogeom vertical diffusion
            nu_min = 1e4
            ne_cgs = ne * M3_TO_CM3
            nHi_cgs = nHi * M3_TO_CM3 if nHi is not None else np.zeros(n_points)
            ne_safe = np.maximum(ne_cgs, 1e-10)
            
            # Initialize sums
            mfp_sum = np.zeros(n_points)
            gr_sum = np.zeros(n_points)
            nu_sum = np.zeros(n_points)
            n_total = np.zeros(n_points)
            
            # H+ contribution
            if nHi is not None and nuHi is not None and THi is not None:
                nuHi_eff = np.maximum(nuHi, nu_min)
                T_eff = Te + THi * nHi_cgs / ne_safe if Te is not None else THi
                mfp_sum += nHi_cgs * 9.79e5 * np.sqrt(np.maximum(T_eff, 0.01) / 1.0) / nuHi_eff
                gr_sum += nHi_cgs * 1.02e2 * np.sqrt(np.maximum(THi, 0.01)) / Br_safe / 1e4
                nu_sum += nHi_cgs * nuHi_eff
                n_total += nHi_cgs
            
            # H2+ contribution
            if nH2i is not None and nuH2i is not None:
                nH2i_cgs = nH2i * M3_TO_CM3
                T_H2i = TH2i if TH2i is not None else (Te if Te is not None else np.full(n_points, 1.0))
                nuH2i_eff = np.maximum(nuH2i, nu_min)
                T_eff = Te + T_H2i * nH2i_cgs / ne_safe if Te is not None else T_H2i
                mfp_sum += nH2i_cgs * 9.79e5 * np.sqrt(np.maximum(T_eff, 0.01) / 2.0) / nuH2i_eff
                gr_sum += nH2i_cgs * 1.02e2 * np.sqrt(2.0 * np.maximum(T_H2i, 0.01)) / Br_safe / 1e4
                nu_sum += nH2i_cgs * nuH2i_eff
                n_total += nH2i_cgs
            
            # H3+ contribution
            if nH3i is not None and nuH3i is not None:
                nH3i_cgs = nH3i * M3_TO_CM3
                T_H3i = TH3i if TH3i is not None else (Te if Te is not None else np.full(n_points, 1.0))
                nuH3i_eff = np.maximum(nuH3i, nu_min)
                T_eff = Te + T_H3i * nH3i_cgs / ne_safe if Te is not None else T_H3i
                mfp_sum += nH3i_cgs * 9.79e5 * np.sqrt(np.maximum(T_eff, 0.01) / 3.0) / nuH3i_eff
                gr_sum += nH3i_cgs * 1.02e2 * np.sqrt(3.0 * np.maximum(T_H3i, 0.01)) / Br_safe / 1e4
                nu_sum += nH3i_cgs * nuH3i_eff
                n_total += nH3i_cgs
            
            # HeII contribution
            if nHeII is not None and nuHeII is not None:
                nHeII_cgs = nHeII * M3_TO_CM3
                T_HeII = THeII if THeII is not None else (Te if Te is not None else np.full(n_points, 1.0))
                nuHeII_eff = np.maximum(nuHeII, nu_min)
                T_eff = Te + T_HeII * nHeII_cgs / ne_safe if Te is not None else T_HeII
                mfp_sum += nHeII_cgs * 9.79e5 * np.sqrt(np.maximum(T_eff, 0.01) / 4.0) / nuHeII_eff
                gr_sum += nHeII_cgs * 1.02e2 * np.sqrt(4.0 * np.maximum(T_HeII, 0.01)) / Br_safe / 1e4
                nu_sum += nHeII_cgs * nuHeII_eff
                n_total += nHeII_cgs
            
            # HeIII contribution
            if nHeIII is not None and nuHeIII is not None:
                nHeIII_cgs = nHeIII * M3_TO_CM3
                T_HeIII = THeIII if THeIII is not None else (Te if Te is not None else np.full(n_points, 1.0))
                nuHeIII_eff = np.maximum(nuHeIII, nu_min)
                T_eff = Te + T_HeIII * nHeIII_cgs / ne_safe if Te is not None else T_HeIII
                mfp_sum += nHeIII_cgs * 9.79e5 * np.sqrt(2.0 * np.maximum(T_eff, 0.01) / 4.0) / nuHeIII_eff
                gr_sum += nHeIII_cgs * 1.02e2 / 2.0 * np.sqrt(4.0 * np.maximum(T_HeIII, 0.01)) / Br_safe / 1e4
                nu_sum += nHeIII_cgs * nuHeIII_eff
                n_total += nHeIII_cgs
            
            # Weighted averages
            n_total_safe = np.maximum(n_total, 1e-10)
            mfp_avg = mfp_sum / n_total_safe
            gr_avg = gr_sum / n_total_safe
            nu_avg = nu_sum / n_total_safe
            
            Dv = Dfsave * 0.333 * nu_avg * mfp_avg * (gr_avg + mfp_avg * Bv / Br_safe)
        else:
            # Fallback
            if D_ion is not None:
                Dv = D_ion * 1e4
            else:
                Dv = np.full(n_points, 1e2)  # 1 m²/s minimum
        
        # Loss time: τ = b² / (2 * Dv)
        Dv_safe = np.maximum(Dv, 1e-10)
        tau_bpol = (b ** 2) / (2.0 * Dv_safe)
        tau_bpol = np.maximum(tau_bpol, 1e-10)
        
        # Add bpol contribution to loss rates
        # For energy, apply gEd factor: τ_E_bpol = τ_bpol / gEd
        k_n += 1.0 / tau_bpol
        k_E += gEd / tau_bpol
    
    return k_n, k_E


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


def compute_limiter_losses(
    R_positions: np.ndarray,
    lHFS: float,
    lLFS: float,
    D_perp: np.ndarray,
    lambda_n: float,
    lambda_E: float,
    ne: np.ndarray,
    Te: np.ndarray,
    nHi: np.ndarray,
    THi: np.ndarray,
    nH2i: np.ndarray = None,
    TH2i: np.ndarray = None,
    nH3i: np.ndarray = None,
    TH3i: np.ndarray = None,
    nHeII: np.ndarray = None,
    THeII: np.ndarray = None,
    nHeIII: np.ndarray = None,
    THeIII: np.ndarray = None,
    Ta0: float = 0.026
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """
    Compute parallel transport losses to limiters using diffusion-based formula.
    
    In the SOL (R < lHFS or R > lLFS), ions are lost along field lines with
    characteristic time based on perpendicular diffusion:
    
        D_⊥ / λ² = 1 / τ_∥   =>   τ = λ² / D_⊥
    
    Separate decay lengths are used for density and energy:
        τ_n = λ_n² / D_⊥  (density loss time)
        τ_E = λ_E² / D_⊥  (energy loss time)
    
    Parameters
    ----------
    R_positions : np.ndarray
        Radial positions [m].
    lHFS : float
        HFS limiter position [m].
    lLFS : float
        LFS limiter position [m].
    D_perp : np.ndarray
        Perpendicular diffusion coefficient [m²/s].
    lambda_n : float
        Density decay length [m].
    lambda_E : float
        Energy decay length [m].
    ne, Te : np.ndarray
        Electron density [m^-3] and temperature [eV].
    nHi, THi : np.ndarray
        H+ density [m^-3] and temperature [eV].
    nH2i, TH2i : np.ndarray, optional
        H2+ density [m^-3] and temperature [eV].
    nH3i, TH3i : np.ndarray, optional
        H3+ density [m^-3] and temperature [eV].
    nHeII, THeII : np.ndarray, optional
        He+ density [m^-3] and temperature [eV].
    nHeIII, THeIII : np.ndarray, optional
        He++ density [m^-3] and temperature [eV].
    Ta0 : float
        Wall/recycled neutral temperature [eV]. Default 0.026 eV (room temp).
        
    Returns
    -------
    dn : Dict[str, np.ndarray]
        Density source terms [m^-3/s].
    dE : Dict[str, np.ndarray]
        Energy source terms [eV·m^-3/s].
    """
    # Initialize output dictionaries
    dn = {}
    dE = {}
    
    # Find SOL mask: points outside the confined region
    # Exclude boundary points (first and last) since they already have decay length BCs
    sol_mask = (R_positions < lHFS) | (R_positions > lLFS)
    sol_mask[0] = False   # HFS boundary - handled by decay length BC
    sol_mask[-1] = False  # LFS boundary - handled by decay length BC
    
    # If no points in SOL, return empty sources
    if not np.any(sol_mask):
        return dn, dE
    
    # Compute loss times from diffusion: τ = λ² / D_⊥
    D_safe = np.maximum(D_perp, 1e-10)  # Prevent division by zero
    tau_n = lambda_n**2 / D_safe  # Density loss time [s]
    tau_E = lambda_E**2 / D_safe  # Energy loss time [s]
    
    # Ensure minimum loss time
    tau_n = np.maximum(tau_n, 1e-10)
    tau_E = np.maximum(tau_E, 1e-10)
    
    # Initialize sources
    dn['e'] = np.zeros_like(ne)
    dE['e'] = np.zeros_like(ne)
    dn['Hi'] = np.zeros_like(ne)
    dE['Hi'] = np.zeros_like(ne)
    
    # --- H+ losses ---
    Z = 1.0
    loss_rate_n_Hi = np.zeros_like(ne)
    loss_rate_E_Hi = np.zeros_like(ne)
    loss_rate_n_Hi[sol_mask] = nHi[sol_mask] / tau_n[sol_mask]
    loss_rate_E_Hi[sol_mask] = (nHi * 1.5 * THi)[sol_mask] / tau_E[sol_mask]
    
    dn['Hi'] -= loss_rate_n_Hi
    dE['Hi'] -= loss_rate_E_Hi
    dn['e'] -= Z * loss_rate_n_Hi
    dE['e'] -= Z * (ne * 1.5 * Te / tau_E)
    dE['e'][~sol_mask] = 0.0
    
    # --- H2+ losses ---
    if nH2i is not None:
        if TH2i is None:
            TH2i = Te.copy()
        
        dn['H2i'] = np.zeros_like(ne)
        dE['H2i'] = np.zeros_like(ne)
        
        Z = 1.0
        loss_rate_n_H2i = np.zeros_like(ne)
        loss_rate_E_H2i = np.zeros_like(ne)
        loss_rate_n_H2i[sol_mask] = nH2i[sol_mask] / tau_n[sol_mask]
        loss_rate_E_H2i[sol_mask] = (nH2i * 1.5 * TH2i)[sol_mask] / tau_E[sol_mask]
        
        dn['H2i'] -= loss_rate_n_H2i
        dE['H2i'] -= loss_rate_E_H2i
        dn['e'] -= Z * loss_rate_n_H2i
        # Electron energy loss already accounted for in H+ section
    
    # --- H3+ losses ---
    if nH3i is not None:
        if TH3i is None:
            TH3i = Te.copy()
        
        dn['H3i'] = np.zeros_like(ne)
        dE['H3i'] = np.zeros_like(ne)
        
        Z = 1.0
        loss_rate_n_H3i = np.zeros_like(ne)
        loss_rate_E_H3i = np.zeros_like(ne)
        loss_rate_n_H3i[sol_mask] = nH3i[sol_mask] / tau_n[sol_mask]
        loss_rate_E_H3i[sol_mask] = (nH3i * 1.5 * TH3i)[sol_mask] / tau_E[sol_mask]
        
        dn['H3i'] -= loss_rate_n_H3i
        dE['H3i'] -= loss_rate_E_H3i
        dn['e'] -= Z * loss_rate_n_H3i
    
    # --- He+ (HeII) losses ---
    if nHeII is not None:
        if THeII is None:
            THeII = Te.copy()
        
        dn['HeII'] = np.zeros_like(ne)
        dE['HeII'] = np.zeros_like(ne)
        
        Z = 1.0
        loss_rate_n_HeII = np.zeros_like(ne)
        loss_rate_E_HeII = np.zeros_like(ne)
        loss_rate_n_HeII[sol_mask] = nHeII[sol_mask] / tau_n[sol_mask]
        loss_rate_E_HeII[sol_mask] = (nHeII * 1.5 * THeII)[sol_mask] / tau_E[sol_mask]
        
        dn['HeII'] -= loss_rate_n_HeII
        dE['HeII'] -= loss_rate_E_HeII
        dn['e'] -= Z * loss_rate_n_HeII
    
    # --- He++ (HeIII) losses ---
    if nHeIII is not None:
        if THeIII is None:
            THeIII = Te.copy()
        
        dn['HeIII'] = np.zeros_like(ne)
        dE['HeIII'] = np.zeros_like(ne)
        
        Z = 2.0
        loss_rate_n_HeIII = np.zeros_like(ne)
        loss_rate_E_HeIII = np.zeros_like(ne)
        loss_rate_n_HeIII[sol_mask] = nHeIII[sol_mask] / tau_n[sol_mask]
        loss_rate_E_HeIII[sol_mask] = (nHeIII * 1.5 * THeIII)[sol_mask] / tau_E[sol_mask]
        
        dn['HeIII'] -= loss_rate_n_HeIII
        dE['HeIII'] -= loss_rate_E_HeIII
        dn['e'] -= Z * loss_rate_n_HeIII
    
    return dn, dE


# def compute_limiter_losses(
#     R_positions: np.ndarray,
#     lHFS: float,
#     lLFS: float,
#     nlimiters: float,
#     ne: np.ndarray,
#     Te: np.ndarray,
#     nHi: np.ndarray,
#     THi: np.ndarray,
#     nH2i: np.ndarray = None,
#     TH2i: np.ndarray = None,
#     nH3i: np.ndarray = None,
#     TH3i: np.ndarray = None,
#     nHeII: np.ndarray = None,
#     THeII: np.ndarray = None,
#     nHeIII: np.ndarray = None,
#     THeIII: np.ndarray = None,
#     Ta0: float = 0.026
# ) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
#     """
#     Compute parallel transport losses to limiters in the SOL region.
    
#     In the SOL (R < lHFS or R > lLFS), ions are lost along field lines to 
#     limiters with characteristic time:
#         tau = 0.66 * 2*pi*R / nlimiters / (cs * ndamp(ne) * ndamp(ni))
#     where cs = 9.79e5 * sqrt(Z/mu * (Te + Ti*ni/ne)) cm/s
    
#     Lost ions recycle as neutrals at wall temperature Ta0.
    
#     Parameters
#     ----------
#     R_positions : np.ndarray
#         Radial positions [m].
#     lHFS : float
#         HFS limiter position [m].
#     lLFS : float
#         LFS limiter position [m].
#     nlimiters : float
#         Number of poloidal limiters.
#     ne, Te : np.ndarray
#         Electron density [m^-3] and temperature [eV].
#     nHi, THi : np.ndarray
#         H+ density [m^-3] and temperature [eV].
#     nH2i, TH2i : np.ndarray, optional
#         H2+ density [m^-3] and temperature [eV].
#     nH3i, TH3i : np.ndarray, optional
#         H3+ density [m^-3] and temperature [eV].
#     nHeII, THeII : np.ndarray, optional
#         He+ density [m^-3] and temperature [eV].
#     nHeIII, THeIII : np.ndarray, optional
#         He++ density [m^-3] and temperature [eV].
#     Ta0 : float
#         Wall/recycled neutral temperature [eV]. Default 0.026 eV (room temp).
        
#     Returns
#     -------
#     dn : Dict[str, np.ndarray]
#         Density source terms [m^-3/s].
#     dE : Dict[str, np.ndarray]
#         Energy source terms [eV·m^-3/s].
#     """
#     # Initialize output dictionaries
#     dn = {}
#     dE = {}
    
#     # Find SOL mask: points outside the confined region
#     # Exclude boundary points (first and last) since they already have decay length BCs
#     sol_mask = (R_positions < lHFS) | (R_positions > lLFS)
#     sol_mask[0] = False   # HFS boundary - handled by decay length BC
#     sol_mask[-1] = False  # LFS boundary - handled by decay length BC
    
#     # If no points in SOL, return empty sources
#     if not np.any(sol_mask):
#         return dn, dE
    
#     # Convert to CGS for rate calculation
#     ne_cgs = ne * M3_TO_CM3
#     nHi_cgs = nHi * M3_TO_CM3
    
#     # R positions in cm for tau calculation
#     R_cm = R_positions * 100.0
    
#     # Energy coefficient for limiters (gElim = 1.0 in C++)
#     gElim = 1.0
    
#     # Initialize sources
#     dn['e'] = np.zeros_like(ne)
#     dE['e'] = np.zeros_like(ne)
#     dn['Hi'] = np.zeros_like(ne)
#     dE['Hi'] = np.zeros_like(ne)
#     dn['H2'] = np.zeros_like(ne)
#     dE['H2'] = np.zeros_like(ne)
    
#     # --- H+ losses ---
#     Z = 1.0
#     mu = 1.0
#     # Sound speed: cs = 9.79e5 * sqrt(Z/mu * (Te + Ti*ni/ne)) [cm/s]
#     T_eff = Te + THi * nHi / np.maximum(ne, 1e-10)
#     cs = 9.79e5 * np.sqrt(Z / mu * np.maximum(T_eff, 0.01))
    
#     # Parallel loss time: tau = 0.66 * 2*pi*R / nlimiters / (cs * ndamp(ne) * ndamp(ni))
#     tau_Hi = np.ones_like(ne) * 1e30  # Large value (no loss) by default
#     damp_factor = cs * ndamp(ne_cgs) * ndamp(nHi_cgs)
#     valid = (damp_factor > 0) & sol_mask
#     tau_Hi[valid] = 0.66 * 2.0 * np.pi * R_cm[valid] / nlimiters / damp_factor[valid]
#     tau_Hi = np.maximum(tau_Hi, 1e-10)  # Prevent division by zero
    
#     # Apply losses only in SOL
#     loss_rate_Hi = nHi / tau_Hi
#     loss_rate_Hi[~sol_mask] = 0.0
    
#     dn['Hi'][sol_mask] -= loss_rate_Hi[sol_mask]
#     dE['Hi'][sol_mask] -= (loss_rate_Hi * 1.5 * THi)[sol_mask]
#     dn['e'][sol_mask] -= Z * loss_rate_Hi[sol_mask]
#     dE['e'][sol_mask] -= (Z * loss_rate_Hi * 1.5 * Te)[sol_mask]
    
#     # # Recycle H+ -> 0.5 H2 at wall temperature
#     # dn['H2'][sol_mask] += 0.5 * loss_rate_Hi[sol_mask]
#     # dE['H2'][sol_mask] += 0.5 * loss_rate_Hi[sol_mask] * 1.5 * Ta0
    
#     # --- H2+ losses ---
#     if nH2i is not None:
#         nH2i_cgs = nH2i * M3_TO_CM3
#         if TH2i is None:
#             TH2i = Te.copy()
        
#         dn['H2i'] = np.zeros_like(ne)
#         dE['H2i'] = np.zeros_like(ne)
        
#         Z = 1.0
#         mu = 2.0
#         T_eff = Te + TH2i * nH2i / np.maximum(ne, 1e-10)
#         cs = 9.79e5 * np.sqrt(Z / mu * np.maximum(T_eff, 0.01))
        
#         tau_H2i = np.ones_like(ne) * 1e30
#         damp_factor = cs * ndamp(ne_cgs) * ndamp(nH2i_cgs)
#         valid = (damp_factor > 0) & sol_mask
#         tau_H2i[valid] = 0.66 * 2.0 * np.pi * R_cm[valid] / nlimiters / damp_factor[valid]
#         tau_H2i = np.maximum(tau_H2i, 1e-10)
        
#         loss_rate_H2i = nH2i / tau_H2i
#         loss_rate_H2i[~sol_mask] = 0.0
        
#         dn['H2i'][sol_mask] -= loss_rate_H2i[sol_mask]
#         dE['H2i'][sol_mask] -= (loss_rate_H2i * 1.5 * TH2i)[sol_mask]
#         dn['e'][sol_mask] -= Z * loss_rate_H2i[sol_mask]
#         dE['e'][sol_mask] -= (Z * loss_rate_H2i * 1.5 * Te)[sol_mask]
        
#         # # Recycle H2+ -> H2 at wall temperature
#         # dn['H2'][sol_mask] += loss_rate_H2i[sol_mask]
#         # dE['H2'][sol_mask] += loss_rate_H2i[sol_mask] * 1.5 * Ta0
    
#     # --- H3+ losses ---
#     if nH3i is not None:
#         nH3i_cgs = nH3i * M3_TO_CM3
#         if TH3i is None:
#             TH3i = Te.copy()
        
#         dn['H3i'] = np.zeros_like(ne)
#         dE['H3i'] = np.zeros_like(ne)
        
#         Z = 1.0
#         mu = 3.0
#         T_eff = Te + TH3i * nH3i / np.maximum(ne, 1e-10)
#         cs = 9.79e5 * np.sqrt(Z / mu * np.maximum(T_eff, 0.01))
        
#         tau_H3i = np.ones_like(ne) * 1e30
#         damp_factor = cs * ndamp(ne_cgs) * ndamp(nH3i_cgs)
#         valid = (damp_factor > 0) & sol_mask
#         tau_H3i[valid] = 0.66 * 2.0 * np.pi * R_cm[valid] / nlimiters / damp_factor[valid]
#         tau_H3i = np.maximum(tau_H3i, 1e-10)
        
#         loss_rate_H3i = nH3i / tau_H3i
#         loss_rate_H3i[~sol_mask] = 0.0
        
#         dn['H3i'][sol_mask] -= loss_rate_H3i[sol_mask]
#         dE['H3i'][sol_mask] -= (loss_rate_H3i * 1.5 * TH3i)[sol_mask]
#         dn['e'][sol_mask] -= Z * loss_rate_H3i[sol_mask]
#         dE['e'][sol_mask] -= (Z * loss_rate_H3i * 1.5 * Te)[sol_mask]
        
#         # # Recycle H3+ -> 1.5 H2 at wall temperature
#         # dn['H2'][sol_mask] += 1.5 * loss_rate_H3i[sol_mask]
#         # dE['H2'][sol_mask] += 1.5 * loss_rate_H3i[sol_mask] * 1.5 * Ta0
    
#     # --- He+ (HeII) losses ---
#     if nHeII is not None:
#         nHeII_cgs = nHeII * M3_TO_CM3
#         if THeII is None:
#             THeII = Te.copy()
        
#         dn['HeII'] = np.zeros_like(ne)
#         dE['HeII'] = np.zeros_like(ne)
#         dn['HeI'] = np.zeros_like(ne)
#         dE['HeI'] = np.zeros_like(ne)
        
#         Z = 1.0
#         mu = 4.0
#         T_eff = Te + THeII * nHeII / np.maximum(ne, 1e-10)
#         cs = 9.79e5 * np.sqrt(Z / mu * np.maximum(T_eff, 0.01))
        
#         tau_HeII = np.ones_like(ne) * 1e30
#         damp_factor = cs * ndamp(ne_cgs) * ndamp(nHeII_cgs)
#         valid = (damp_factor > 0) & sol_mask
#         tau_HeII[valid] = 0.66 * 2.0 * np.pi * R_cm[valid] / nlimiters / damp_factor[valid]
#         tau_HeII = np.maximum(tau_HeII, 1e-10)
        
#         loss_rate_HeII = nHeII / tau_HeII
#         loss_rate_HeII[~sol_mask] = 0.0
        
#         dn['HeII'][sol_mask] -= loss_rate_HeII[sol_mask]
#         dE['HeII'][sol_mask] -= (loss_rate_HeII * 1.5 * THeII)[sol_mask]
#         dn['e'][sol_mask] -= Z * loss_rate_HeII[sol_mask]
#         dE['e'][sol_mask] -= (Z * loss_rate_HeII * 1.5 * Te)[sol_mask]
        
#         # # Recycle He+ -> HeI at wall temperature
#         # dn['HeI'][sol_mask] += loss_rate_HeII[sol_mask]
#         # dE['HeI'][sol_mask] += loss_rate_HeII[sol_mask] * 1.5 * Ta0
    
#     # --- He++ (HeIII) losses ---
#     if nHeIII is not None:
#         nHeIII_cgs = nHeIII * M3_TO_CM3
#         if THeIII is None:
#             THeIII = Te.copy()
        
#         if 'HeIII' not in dn:
#             dn['HeIII'] = np.zeros_like(ne)
#             dE['HeIII'] = np.zeros_like(ne)
#         if 'HeI' not in dn:
#             dn['HeI'] = np.zeros_like(ne)
#             dE['HeI'] = np.zeros_like(ne)
        
#         Z = 2.0
#         mu = 4.0
#         T_eff = Te + THeIII * nHeIII / np.maximum(ne, 1e-10)
#         cs = 9.79e5 * np.sqrt(Z / mu * np.maximum(T_eff, 0.01))
        
#         tau_HeIII = np.ones_like(ne) * 1e30
#         damp_factor = cs * ndamp(ne_cgs) * ndamp(nHeIII_cgs)
#         valid = (damp_factor > 0) & sol_mask
#         tau_HeIII[valid] = 0.66 * 2.0 * np.pi * R_cm[valid] / nlimiters / damp_factor[valid]
#         tau_HeIII = np.maximum(tau_HeIII, 1e-10)
        
#         loss_rate_HeIII = nHeIII / tau_HeIII
#         loss_rate_HeIII[~sol_mask] = 0.0
        
#         dn['HeIII'][sol_mask] -= loss_rate_HeIII[sol_mask]
#         dE['HeIII'][sol_mask] -= (loss_rate_HeIII * 1.5 * THeIII)[sol_mask]
#         dn['e'][sol_mask] -= Z * loss_rate_HeIII[sol_mask]
#         dE['e'][sol_mask] -= (Z * loss_rate_HeIII * 1.5 * Te)[sol_mask]
        
#         # # Recycle He++ -> HeI at wall temperature
#         # dn['HeI'][sol_mask] += loss_rate_HeIII[sol_mask]
#         # dE['HeI'][sol_mask] += loss_rate_HeIII[sol_mask] * 1.5 * Ta0
    
#     return dn, dE

def compute_bpol_losses(
    Br: np.ndarray,
    Bv: float,
    b: float,
    ne: np.ndarray,
    Te: np.ndarray,
    nHi: np.ndarray,
    THi: np.ndarray,
    EHi: np.ndarray,
    nuHi: np.ndarray,
    nH2i: np.ndarray = None,
    TH2i: np.ndarray = None,
    EH2i: np.ndarray = None,
    nuH2i: np.ndarray = None,
    nH3i: np.ndarray = None,
    TH3i: np.ndarray = None,
    EH3i: np.ndarray = None,
    nuH3i: np.ndarray = None,
    nHeII: np.ndarray = None,
    THeII: np.ndarray = None,
    EHeII: np.ndarray = None,
    nuHeII: np.ndarray = None,
    nHeIII: np.ndarray = None,
    THeIII: np.ndarray = None,
    EHeIII: np.ndarray = None,
    nuHeIII: np.ndarray = None,
    Dfsave: float = 1.0,
    D_ion: np.ndarray = None,
    diffusion_model: str = 'gyrogeom',
    gEd: float = 1.0,
    Ta0: float = 0.026
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """
    Compute vertical diffusion losses along inclined field lines (bpol_function).
    
    This implements vertical transport losses where particles diffuse along 
    the poloidal magnetic field direction. The diffusion coefficient accounts
    for cross-field transport enhanced by the vertical field component.
    
    From C++ Tomator1D bpol_function (lines 671-844):
        - Computes density-weighted mean free path and gyroradius for all ions
        - Vertical diffusion: Dv = Dfsave * 0.333 * nu * mfp * (gr + mfp * Bv/Br)
        - Loss time: tau = b² / (2 * Dv)
    
    Parameters
    ----------
    Br : np.ndarray
        Local toroidal magnetic field [T]. Varies with position as Bt * R0 / R.
    Bv : float
        Vertical magnetic field component [T].
    b : float
        Vertical extent of plasma [cm].
    ne, Te : np.ndarray
        Electron density [m^-3] and temperature [eV].
    nHi, THi, EHi : np.ndarray
        H+ density [m^-3], temperature [eV], and energy density [eV·m^-3].
    nuHi : np.ndarray
        H+ collision frequency [s^-1].
    nH2i, TH2i, EH2i, nuH2i : np.ndarray, optional
        H2+ quantities.
    nH3i, TH3i, EH3i, nuH3i : np.ndarray, optional
        H3+ quantities.
    nHeII, THeII, EHeII, nuHeII : np.ndarray, optional
        He+ quantities.
    nHeIII, THeIII, EHeIII, nuHeIII : np.ndarray, optional
        He++ quantities.
    Dfsave : float
        Diffusion scaling factor. Default 1.0.
    D_ion : np.ndarray, optional
        Pre-computed ion diffusion coefficient [m²/s]. Used directly for
        'fixed' and 'bohm' models. For 'gyrogeom', this is ignored and
        vertical diffusion is computed with Bv/Br formula.
    diffusion_model : str
        Diffusion model type: 'fixed', 'bohm', or 'gyrogeom'. Default 'gyrogeom'.
        - 'fixed' or 'bohm': Use D_ion directly (no Bv dependence)
        - 'gyrogeom': Compute Dv = Dfsave * 0.333 * nu * mfp * (gr + mfp * Bv/Br)
    gEd : float
        Energy loss rate factor. Default 1.0.
    Ta0 : float
        Wall temperature for recycled neutrals [eV]. Default 0.026.
        
    Returns
    -------
    dn : Dict[str, np.ndarray]
        Density source terms [m^-3/s].
    dE : Dict[str, np.ndarray]
        Energy source terms [eV·m^-3/s].
    """
    # Initialize output dictionaries
    dn = {}
    dE = {}
    
    # Convert densities to CGS (cm^-3)
    ne_cgs = ne * M3_TO_CM3
    nHi_cgs = nHi * M3_TO_CM3
    
    # Initialize collision frequencies and source arrays
    n_points = len(ne)
    
    dn['e'] = np.zeros(n_points)
    dE['e'] = np.zeros(n_points)
    dn['Hi'] = np.zeros(n_points)
    dE['Hi'] = np.zeros(n_points)
    
    # Minimum collision frequency to prevent numerical issues
    # C++ uses max(nu, 1e4) in transport.cpp line 72
    nu_min = 1e4  # [s^-1]
    
    # Br is in Tesla, formulas use Gauss (1 T = 1e4 Gauss)
    # In C++: gr = ... / Br[im] / 1e4, equivalent to dividing by B in Gauss
    Br_safe = np.maximum(np.abs(Br), 1e-6)  # Prevent division by zero [T]
    
    # Compute vertical diffusion coefficient based on diffusion model
    # - Fixed/Bohm: use pre-computed D_ion directly (already calculated for radial transport)
    # - Gyrogeom: compute vertical-specific formula with Bv/Br term
    
    if diffusion_model.lower() in ('fixed', 'bohm') and D_ion is not None:
        # Use pre-computed diffusion coefficient directly
        # Convert from m²/s to cm²/s (factor 1e4) since b is in cm
        Dv = D_ion * 1e4  # [cm²/s]
    elif diffusion_model.lower() == 'gyrogeom':
        # Compute weighted mean free path, gyroradius, and collision frequency
        # mfp = n * 9.79e5 * sqrt((Te + Ti*ni/ne) / mu) / max(nu, 1e4) [cm]
        # gr = n * 1.02e2 / Z * sqrt(mu * Ti) / Br / 1e4 [cm]
        # nu_weighted = n * max(nu, 1e4)
        
        ne_safe = np.maximum(ne_cgs, 1e-10)
        
        # H+ contribution (mu=1, Z=1)
        nuHi_eff = np.maximum(nuHi, nu_min)
        T_eff_Hi = Te + THi * nHi_cgs / ne_safe
        mfp_sum = nHi_cgs * 9.79e5 * np.sqrt(np.maximum(T_eff_Hi, 0.01) / 1.0) / nuHi_eff
        gr_sum = nHi_cgs * 1.02e2 / 1.0 * np.sqrt(1.0 * np.maximum(THi, 0.01)) / Br_safe / 1e4
        nu_sum = nHi_cgs * nuHi_eff
        n_total = nHi_cgs.copy()
        
        # H2+ contribution (mu=2, Z=1)
        if nH2i is not None and nuH2i is not None:
            nH2i_cgs = nH2i * M3_TO_CM3
            if TH2i is None:
                TH2i = Te.copy()
            nuH2i_eff = np.maximum(nuH2i, nu_min)
            T_eff_H2i = Te + TH2i * nH2i_cgs / ne_safe
            mfp_sum += nH2i_cgs * 9.79e5 * np.sqrt(np.maximum(T_eff_H2i, 0.01) / 2.0) / nuH2i_eff
            gr_sum += nH2i_cgs * 1.02e2 / 1.0 * np.sqrt(2.0 * np.maximum(TH2i, 0.01)) / Br_safe / 1e4
            nu_sum += nH2i_cgs * nuH2i_eff
            n_total += nH2i_cgs
        
        # H3+ contribution (mu=3, Z=1)
        if nH3i is not None and nuH3i is not None:
            nH3i_cgs = nH3i * M3_TO_CM3
            if TH3i is None:
                TH3i = Te.copy()
            nuH3i_eff = np.maximum(nuH3i, nu_min)
            T_eff_H3i = Te + TH3i * nH3i_cgs / ne_safe
            mfp_sum += nH3i_cgs * 9.79e5 * np.sqrt(np.maximum(T_eff_H3i, 0.01) / 3.0) / nuH3i_eff
            gr_sum += nH3i_cgs * 1.02e2 / 1.0 * np.sqrt(3.0 * np.maximum(TH3i, 0.01)) / Br_safe / 1e4
            nu_sum += nH3i_cgs * nuH3i_eff
            n_total += nH3i_cgs
        
        # HeII contribution (mu=4, Z=1)
        if nHeII is not None and nuHeII is not None:
            nHeII_cgs = nHeII * M3_TO_CM3
            if THeII is None:
                THeII = Te.copy()
            nuHeII_eff = np.maximum(nuHeII, nu_min)
            T_eff_HeII = Te + THeII * nHeII_cgs / ne_safe
            mfp_sum += nHeII_cgs * 9.79e5 * np.sqrt(np.maximum(T_eff_HeII, 0.01) / 4.0) / nuHeII_eff
            gr_sum += nHeII_cgs * 1.02e2 / 1.0 * np.sqrt(4.0 * np.maximum(THeII, 0.01)) / Br_safe / 1e4
            nu_sum += nHeII_cgs * nuHeII_eff
            n_total += nHeII_cgs
        
        # HeIII contribution (mu=4, Z=2)
        if nHeIII is not None and nuHeIII is not None:
            nHeIII_cgs = nHeIII * M3_TO_CM3
            if THeIII is None:
                THeIII = Te.copy()
            nuHeIII_eff = np.maximum(nuHeIII, nu_min)
            # Note: Z=2 appears in sqrt factor for mfp, and /Z for gyroradius
            T_eff_HeIII = Te + THeIII * nHeIII_cgs / ne_safe
            mfp_sum += nHeIII_cgs * 9.79e5 * np.sqrt(2.0 * np.maximum(T_eff_HeIII, 0.01) / 4.0) / nuHeIII_eff
            gr_sum += nHeIII_cgs * 1.02e2 / 2.0 * np.sqrt(4.0 * np.maximum(THeIII, 0.01)) / Br_safe / 1e4
            nu_sum += nHeIII_cgs * nuHeIII_eff
            n_total += nHeIII_cgs
        
        # Compute density-weighted averages
        n_total_safe = np.maximum(n_total, 1e-10)
        mfp_avg = mfp_sum / n_total_safe  # [cm]
        gr_avg = gr_sum / n_total_safe    # [cm]
        nu_avg = nu_sum / n_total_safe    # [s^-1]
        
        # Vertical diffusion coefficient: Dv = Dfsave * 0.333 * nu * mfp * (gr + mfp * Bv/Br) [cm^2/s]
        Dv = Dfsave * 0.333 * nu_avg * mfp_avg * (gr_avg + mfp_avg * Bv / Br_safe)
    else:
        # Fallback: if D_ion provided without model specification, use it
        if D_ion is not None:
            Dv = D_ion * 1e4  # [cm²/s]
        else:
            # Default to small diffusion if nothing provided
            Dv = np.full(n_points, 1e2)  # 1e2 cm²/s = 0.01 m²/s minimum
    
    # Loss time: tau = b^2 / (2 * Dv) [s]
    Dv_safe = np.maximum(Dv, 1e-10)
    tau_ion = (b ** 2) / (2.0 * Dv_safe)
    tau_ion = np.maximum(tau_ion, 1e-10)  # Prevent division by zero
    
    # --- H+ losses ---
    Z = 1.0
    dn['Hi'] -= nHi / tau_ion
    dE['Hi'] -= EHi / (tau_ion / gEd)
    dn['e'] -= Z * nHi / tau_ion
    dE['e'] -= Z * nHi / (tau_ion / gEd) * Te * 1.5
    
    # --- H2+ losses ---
    if nH2i is not None:
        dn['H2i'] = np.zeros(n_points)
        dE['H2i'] = np.zeros(n_points)
        if EH2i is None:
            EH2i = 1.5 * TH2i * nH2i
        
        Z = 1.0
        dn['H2i'] -= nH2i / tau_ion
        dE['H2i'] -= EH2i / (tau_ion / gEd)
        dn['e'] -= Z * nH2i / tau_ion
        dE['e'] -= Z * nH2i / (tau_ion / gEd) * Te * 1.5
    
    # --- H3+ losses ---
    if nH3i is not None:
        dn['H3i'] = np.zeros(n_points)
        dE['H3i'] = np.zeros(n_points)
        if EH3i is None:
            EH3i = 1.5 * TH3i * nH3i
        
        Z = 1.0
        dn['H3i'] -= nH3i / tau_ion
        dE['H3i'] -= EH3i / (tau_ion / gEd)
        dn['e'] -= Z * nH3i / tau_ion
        dE['e'] -= Z * nH3i / (tau_ion / gEd) * Te * 1.5
    
    # --- HeII (He+) losses ---
    if nHeII is not None:
        dn['HeII'] = np.zeros(n_points)
        dE['HeII'] = np.zeros(n_points)
        if EHeII is None:
            EHeII = 1.5 * THeII * nHeII
        
        Z = 1.0
        dn['HeII'] -= nHeII / tau_ion
        dE['HeII'] -= EHeII / (tau_ion / gEd)
        dn['e'] -= Z * nHeII / tau_ion
        dE['e'] -= Z * nHeII / (tau_ion / gEd) * Te * 1.5
    
    # --- HeIII (He++) losses ---
    if nHeIII is not None:
        dn['HeIII'] = np.zeros(n_points)
        dE['HeIII'] = np.zeros(n_points)
        if EHeIII is None:
            EHeIII = 1.5 * THeIII * nHeIII
        
        Z = 2.0
        dn['HeIII'] -= nHeIII / tau_ion
        dE['HeIII'] -= EHeIII / (tau_ion / gEd)
        dn['e'] -= Z * nHeIII / tau_ion
        dE['e'] -= Z * nHeIII / (tau_ion / gEd) * Te * 1.5
    
    return dn, dE
