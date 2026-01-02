"""
Coupled power module for RF power deposition.

This module provides functions to compute RF power deposition profiles
for different coupling modes. Power is deposited to electrons and ions
to allow steady-state plasma profiles to develop.

Ported from C++ Tomator1D coupledpower.cpp and related files.
"""

import numpy as np
from typing import Dict, Any, Optional, Tuple
from dataclasses import dataclass, field


@dataclass
class PowerDepositionResult:
    """Container for power deposition arrays.
    
    All arrays are in units of [eV/cm³/s] for direct use in energy equations.
    """
    PRFe: np.ndarray      # Power to electrons
    PRFHi: np.ndarray     # Power to H+
    PRFH2i: np.ndarray    # Power to H2+
    PRFH3i: np.ndarray    # Power to H3+
    PRFHeII: np.ndarray   # Power to He+
    PRFHeIII: np.ndarray  # Power to He++
    pecabs: float = 0.0   # Absorbed power fraction (for diagnostics)


@dataclass
class CoupledPowerState:
    """State variables for coupled power controller.
    
    Used to maintain PID controller state between time steps.
    """
    PIerrorP: float = 0.0  # Integral error for PID controller
    pecabs: float = 0.1    # Current absorbed power fraction
    prev_error: float = 0.0  # Previous error for derivative term
    prev_nefact: float = 1.0  # Previous nefact for derivative calculation


class CoupledPower:
    """
    Coupled power manager for RF power deposition.
    
    This class handles the selection and computation of RF power coupling
    to the plasma. Different modes are available:
    
    - nopower: No RF power deposition
    - fixpowerfrac: Fixed power fraction absorbed
    - nefix: PI control to maintain fixed electron density
    - proptone: Power proportional to electron density
    
    Parameters
    ----------
    params : dict
        Simulation parameters containing:
        - Prf : RF power [kW]
        - Rdep : Deposition radius [m]
        - widthech : EC deposition width [m]
        - pecabs0 : Initial absorbed power fraction
        - echbackground : Background power fraction
        - necfix : Target electron density for nefix mode [m^-3]
        - ic : Radial index for nefix density control
        - tauP : PI controller time constant [s]
        - b : Vertical extent [m]
    
    Attributes
    ----------
    mode : str
        Active power coupling mode
    state : CoupledPowerState
        Controller state (PI error, absorbed fraction)
    """
    
    # Available power modes
    MODES = ['nopower', 'fixpowerfrac', 'nefix', 'proptone']
    
    def __init__(self, params: Dict[str, Any]):
        """Initialize coupled power manager.
        
        Parameters
        ----------
        params : dict
            Simulation parameters
        """
        self.params = params
        self.state = CoupledPowerState()
        
        # Parse RF power parameters
        self.Prf = params.get('Prf', 0.0)  # kW
        self.Rdep = params.get('Rdep', 0.9)  # m (convert from cm if needed)
        self.widthech = params.get('widthech', 0.02)  # m
        self.pecabs0 = params.get('pecabs0', 0.1)
        self.echbackground = params.get('echbackground', 1e-5)
        self.b = params.get('b', 0.75)  # m
        
        # nefix mode parameters
        self.necfix = params.get('necfix', 1e19)  # m^-3
        self.ic = params.get('ic', 90)  # control index
        self.tauP = params.get('tauP', 5e-5)  # s
        self.Pini = params.get('Pini', 0.0)  # Initial PIerrorP value
        
        # Initialize state
        self.state.pecabs = self.pecabs0
        self.state.PIerrorP = self.Pini  # Initialize from Pini (C++ behavior)
        
        # Determine mode from boolean flags
        self.mode = self._determine_mode()
    
    def _determine_mode(self) -> str:
        """Determine power coupling mode from simulation flags."""
        if self.params.get('bnopower', False):
            return 'nopower'
        elif self.params.get('bnefix', False):
            return 'nefix'
        elif self.params.get('bfixpowerfrac', False):
            return 'fixpowerfrac'
        elif self.params.get('bproptone', False):
            return 'proptone'
        else:
            # Default to nopower if no mode specified
            return 'nopower'
    
    def set_mode(self, mode: str):
        """Set the power coupling mode.
        
        Parameters
        ----------
        mode : str
            One of 'nopower', 'fixpowerfrac', 'nefix', 'proptone'
        """
        if mode not in self.MODES:
            raise ValueError(f"Unknown power mode: {mode}. Valid modes: {self.MODES}")
        self.mode = mode
    
    def compute_power(
        self,
        R: np.ndarray,
        ne: np.ndarray,
        Te: np.ndarray,
        dt: float,
        t: float,
        nue: Optional[np.ndarray] = None
    ) -> PowerDepositionResult:
        """
        Compute RF power deposition for all species.
        
        Parameters
        ----------
        R : np.ndarray
            Radial coordinates [m]
        ne : np.ndarray
            Electron density [m^-3]
        Te : np.ndarray
            Electron temperature [eV]
        dt : float
            Time step [s]
        t : float
            Current simulation time [s]
        nue : np.ndarray, optional
            Electron collision frequency (for proptone mode) [s^-1]
        
        Returns
        -------
        PowerDepositionResult
            Power deposition arrays for all species [eV/cm³/s]
        """
        nmesh = len(R)
        
        # Select mode
        if self.mode == 'nopower':
            result = self._nopower(nmesh)
        elif self.mode == 'fixpowerfrac':
            result = self._fixpowerfrac(R, nmesh)
        elif self.mode == 'nefix':
            result = self._nefix(R, ne, Te, dt, nmesh)
        elif self.mode == 'proptone':
            result = self._proptone(R, ne, nue, nmesh)
        else:
            result = self._nopower(nmesh)
        
        # Apply power ramp-up at start of simulation
        dtpramp = self.params.get('dtpramp', 0.0)
        if t < 1.02 * dtpramp and dtpramp > 0:
            ramp_factor = (t + 0.02 * dtpramp) / (1.02 * dtpramp)
            result.PRFe *= ramp_factor
            result.PRFHi *= ramp_factor
            result.PRFH2i *= ramp_factor
            result.PRFH3i *= ramp_factor
            result.PRFHeII *= ramp_factor
            result.PRFHeIII *= ramp_factor
        
        return result
    
    def _nopower(self, nmesh: int) -> PowerDepositionResult:
        """No power deposition mode."""
        return PowerDepositionResult(
            PRFe=np.zeros(nmesh),
            PRFHi=np.zeros(nmesh),
            PRFH2i=np.zeros(nmesh),
            PRFH3i=np.zeros(nmesh),
            PRFHeII=np.zeros(nmesh),
            PRFHeIII=np.zeros(nmesh),
            pecabs=0.0
        )
    
    def _compute_ec_profile(
        self, 
        R: np.ndarray,
        ne: Optional[np.ndarray] = None,
        include_background: bool = True
    ) -> np.ndarray:
        """
        Compute normalized EC power deposition profile.
        
        Parameters
        ----------
        R : np.ndarray
            Radial coordinates [m]
        ne : np.ndarray, optional
            Electron density for background profile [m^-3]
        include_background : bool
            Whether to include background profile weighted by ne
        
        Returns
        -------
        np.ndarray
            Normalized power deposition profile
        """
        nmesh = len(R)
        
        # Gaussian EC deposition profile: R * exp(-((R - Rdep)/widthech)^2)
        eprof = R * np.exp(-((R - self.Rdep) / self.widthech)**2)
        
        if include_background and ne is not None:
            # Add background proportional to ne (using minor radius 'a' for width, as in C++)
            a = self.params.get('a', R.max() - R.min())  # minor radius
            eprof_bg = R * np.exp(-((R - self.Rdep) / a)**2)
            eprof += self.echbackground * eprof_bg
        
        # Volume-integrate for normalization: 2*b*pi * integral(eprof * 2*R * dR)
        # Using cylindrical shell volumes: 2*pi*b * (R[i+1]^2 - R[i]^2)
        inteprof = 0.0
        for i in range(nmesh - 1):
            inteprof += (2.0 * self.b * np.pi * 0.5 * 
                        (eprof[i] + eprof[i+1]) * 
                        (R[i+1]**2 - R[i]**2))
        
        # Normalize
        if inteprof > 0:
            eprof = eprof / inteprof
        
        return eprof
    
    def _fixpowerfrac(self, R: np.ndarray, nmesh: int) -> PowerDepositionResult:
        """
        Fixed power fraction mode.
        
        A constant fraction of launched power is absorbed with a 
        Gaussian deposition profile centered at Rdep.
        """
        pecabs = self.pecabs0
        self.state.pecabs = pecabs
        self.state.PIerrorP = 0.0
        
        # Compute normalized EC profile
        eprof = self._compute_ec_profile(R, include_background=True)
        
        # Power deposition: pecabs * 6.24e18 [eV/J] * eprof * (Prf [kW] * 1e3 [W/kW])
        # Result in eV/cm³/s
        power_factor = pecabs * 6.24e18 * (self.Prf * 1e3)
        PRFe = power_factor * eprof
        
        return PowerDepositionResult(
            PRFe=PRFe,
            PRFHi=np.zeros(nmesh),
            PRFH2i=np.zeros(nmesh),
            PRFH3i=np.zeros(nmesh),
            PRFHeII=np.zeros(nmesh),
            PRFHeIII=np.zeros(nmesh),
            pecabs=pecabs
        )
    
    def _nefix(
        self, 
        R: np.ndarray, 
        ne: np.ndarray, 
        Te: np.ndarray,
        dt: float,
        nmesh: int
    ) -> PowerDepositionResult:
        """
        Fixed density mode with PID controller.
        
        Adjusts absorbed power fraction to maintain target electron density
        at the control location (ic). 
        
        PID Controller:
        - P (Proportional): 1/nefact scaling provides immediate response
        - I (Integral): PIerrorP accumulates error to eliminate steady-state offset
        - D (Derivative): Dampens oscillations by opposing rapid changes
        """
        # Convert ne to CGS for consistency with C++ (m^-3 -> cm^-3)
        ne_cgs = ne * 1e-6  # m^-3 -> cm^-3
        necfix_cgs = self.necfix * 1e-6
        
        # Compute density ratio at control location
        ic = min(self.ic, nmesh - 1)
        nefact = ne_cgs[ic] / necfix_cgs if necfix_cgs > 0 else 1.0
        nefact = max(nefact, 0.01)  # Avoid division by zero
        
        # Compute line-integrated electron density for background profile
        inteneprof = 0.0
        for i in range(nmesh - 1):
            inteneprof += (2.0 * self.b * np.pi * 0.5 * 
                         (ne_cgs[i] + ne_cgs[i+1]) * 
                         (R[i+1]**2 - R[i]**2))
        
        # Compute EC deposition profile
        eprof = R * np.exp(-((R - self.Rdep) / self.widthech)**2)
        
        # Integrate EC profile
        inteprof = 0.0
        for i in range(nmesh - 1):
            inteprof += (2.0 * self.b * np.pi * 0.5 * 
                        (eprof[i] + eprof[i+1]) * 
                        (R[i+1]**2 - R[i]**2))
        
        # Combine EC profile with background (ne-weighted)
        if inteprof > 0 and inteneprof > 0:
            eprof = ((1.0 - self.echbackground) * eprof / inteprof + 
                     self.echbackground * ne_cgs / inteneprof)
        
        # Find Te at deposition location
        i = 0
        while i < nmesh - 1 and R[i] < self.Rdep:
            i += 1
        Te_dep = max(Te[i], Te[min(i+1, nmesh-1)])
        Te_dep = min(Te_dep, 5000.0)  # Cap at 5000 eV
        
        # Temperature-dependent absorbed fraction (base value)
        # Log-linear fit: log10(pecabs) = slope * log10(Te) - const
        dTe = np.log10(54.8 / 50.0)
        dpabs = np.log10(0.109052 / 0.100026)
        cst = np.log10(50.0) * dpabs / dTe - np.log10(0.100026)
        pecabs_base = 10.0**(np.log10(max(Te_dep, 1.0)) * dpabs / dTe - cst)
        
        # ===== PID Controller =====
        #
        # Error: e = 1 - nefact
        #   - Positive when density below target (need more power)
        #   - Negative when density above target (need less power)
        
        error = 1.0 - nefact
        
        # PID gains (tuned for stability)
        Kp = 0.5   # Proportional gain
        Ki = 0.3   # Integral gain
        Kd = 2.0   # Derivative gain (on filtered derivative)
        
        # Exponential relaxation factors
        alpha = 1.0 - np.exp(-dt / self.tauP)            # For integral term
        alphaD = 1.0 - np.exp(-dt / (10.0 * self.tauP))  # For derivative filter (slower)
        
        # --- Proportional term ---
        # Use softer scaling: nefact^(-Kp) instead of 1/nefact
        P_term = nefact ** (-Kp)
        
        # --- Integral term ---
        # Accumulate error with relaxation
        self.state.PIerrorP += alpha * Ki * error
        
        # Anti-windup: limit integral term
        PIerrorP_max = 1.5
        PIerrorP_min = -0.6
        self.state.PIerrorP = np.clip(self.state.PIerrorP, PIerrorP_min, PIerrorP_max)
        
        I_term = self.state.PIerrorP
        
        # --- Derivative term ---
        # Use FILTERED derivative to avoid amplifying rapid changes
        # Instead of d(nefact)/dt, use exponentially smoothed change
        # This prevents the derivative from blowing up with small dt
        #
        # Filtered derivative: D_filtered = alphaD * (nefact - prev_nefact) / tauD
        # This gives a smooth derivative that doesn't depend on dt directly
        
        delta_nefact = nefact - self.state.prev_nefact
        
        # Filtered derivative (bounded by tauD, not dt)
        D_term = -Kd * alphaD * delta_nefact
        
        # Store current values for next iteration
        self.state.prev_error = error
        self.state.prev_nefact = nefact
        
        # Combine PID terms
        # P_term provides base scaling
        # I_term eliminates steady-state error  
        # D_term dampens oscillations
        correction = P_term * (1.0 + I_term + D_term)
        correction = np.clip(correction, 0.05, 10.0)  # Reasonable bounds
        
        pecabs = pecabs_base * correction
        
        # Clamp pecabs to valid range [0.001, 1.0]
        pecabs = np.clip(pecabs, 0.001, 1.0)
        
        self.state.pecabs = pecabs
        
        # Power deposition
        power_factor = pecabs * 6.24e18 * (self.Prf * 1e3)
        PRFe = power_factor * eprof
        
        return PowerDepositionResult(
            PRFe=PRFe,
            PRFHi=np.zeros(nmesh),
            PRFH2i=np.zeros(nmesh),
            PRFH3i=np.zeros(nmesh),
            PRFHeII=np.zeros(nmesh),
            PRFHeIII=np.zeros(nmesh),
            pecabs=pecabs
        )
    
    def _proptone(
        self, 
        R: np.ndarray, 
        ne: np.ndarray,
        nue: Optional[np.ndarray],
        nmesh: int
    ) -> PowerDepositionResult:
        """
        Power proportional to electron density mode.
        
        Deposits power proportionally to ne * nue * R (density-weighted
        collision profile). Used for low-density startup phases.
        """
        # Convert ne to CGS
        ne_cgs = ne * 1e-6  # m^-3 -> cm^-3
        
        if nue is None:
            # If no collision frequency provided, use ne only
            nue = np.ones_like(ne)
        
        # Get plasma volume
        Vpl = self.params.get('Vpl', 1.0)  # m³
        Vpl_cgs = Vpl * 1e6  # cm³
        
        # Line-integrated density-weighted collision frequency
        lnenue = 0.0
        for i in range(nmesh - 1):
            lnenue += 0.5 * ((ne_cgs[i] * nue[i] * R[i] + 
                             ne_cgs[i+1] * nue[i+1] * R[i+1]) * 
                            (R[i+1] - R[i]))
        lnenue = lnenue / (R[-1] - R[0])
        
        if lnenue <= 0:
            return self._nopower(nmesh)
        
        # Power deposition proportional to ne * nue * R
        PRFe = (6.24e18 * self.Prf / Vpl_cgs * 
                ne_cgs * nue * R / lnenue)
        
        return PowerDepositionResult(
            PRFe=PRFe,
            PRFHi=np.zeros(nmesh),
            PRFH2i=np.zeros(nmesh),
            PRFH3i=np.zeros(nmesh),
            PRFHeII=np.zeros(nmesh),
            PRFHeIII=np.zeros(nmesh),
            pecabs=self.state.pecabs
        )
    
    def reset_controller(self):
        """Reset PI controller state."""
        self.state.PIerrorP = self.Pini
        self.state.pecabs = self.pecabs0


def compute_coupled_power(
    params: Dict[str, Any],
    R: np.ndarray,
    ne: np.ndarray,
    Te: np.ndarray,
    dt: float,
    t: float,
    mode: Optional[str] = None,
    nue: Optional[np.ndarray] = None,
    state: Optional[CoupledPowerState] = None
) -> Tuple[PowerDepositionResult, CoupledPowerState]:
    """
    Convenience function to compute coupled power in a single call.
    
    Parameters
    ----------
    params : dict
        Simulation parameters
    R : np.ndarray
        Radial coordinates [m]
    ne : np.ndarray
        Electron density [m^-3]
    Te : np.ndarray
        Electron temperature [eV]
    dt : float
        Time step [s]
    t : float
        Current simulation time [s]
    mode : str, optional
        Power mode override (uses params flags if not specified)
    nue : np.ndarray, optional
        Electron collision frequency [s^-1]
    state : CoupledPowerState, optional
        Previous controller state (for continuity)
    
    Returns
    -------
    PowerDepositionResult
        Power deposition arrays
    CoupledPowerState
        Updated controller state
    """
    cp = CoupledPower(params)
    
    # Restore state if provided
    if state is not None:
        cp.state = state
    
    # Override mode if specified
    if mode is not None:
        cp.set_mode(mode)
    
    result = cp.compute_power(R, ne, Te, dt, t, nue)
    
    return result, cp.state
