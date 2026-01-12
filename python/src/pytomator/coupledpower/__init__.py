"""
Coupled power deposition subpackage.

This module implements RF power coupling and deposition calculations,
ported from the C++ Tomator1D coupledpower.cpp.

Available modes:
- nopower: No RF power deposition
- fixpowerfrac: Fixed fraction of launched power absorbed
- nefix: PI controller to maintain fixed electron density
- proptone: Power proportional to electron density
"""

from .coupledpower import (
    CoupledPower,
    CoupledPowerState,
    PowerDepositionResult,
    compute_coupled_power,
)

__all__ = [
    'CoupledPower',
    'CoupledPowerState', 
    'PowerDepositionResult',
    'compute_coupled_power',
]

__all__ = []
