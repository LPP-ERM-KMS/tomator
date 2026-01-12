"""
HYDHEL database parser for reaction rate polynomial coefficients.

This module reads the HYDHEL.tex file and extracts polynomial fit coefficients
for various plasma reactions. The coefficients are used to calculate rate
coefficients as:

    k(T) = exp(sum_i b_i * (ln T)^i)  for T in eV

The file is parsed once at module import time and cached.
"""

import os
import re
import numpy as np
from typing import Dict, Optional, Tuple
from pathlib import Path


# Default path to HYDHEL file (relative to this file's location)
# hydhel.py is at: tomator/python/tomator_dolfinx/reactions/hydhel.py
# hydhel.tex is at: tomator/src/SimParams/Public/hydhel.tex
# So we go: reactions/ -> tomator_dolfinx/ -> python/ -> tomator/ -> src/
_DEFAULT_HYDHEL_PATH = Path(__file__).parent.parent.parent.parent / "src" / "SimParams" / "Public" / "hydhel.tex"

# Cache for loaded coefficients
_HYDHEL_COEFFICIENTS: Dict[str, np.ndarray] = {}
_HYDHEL_METADATA: Dict[str, Dict] = {}
_HYDHEL_LOADED = False


def _parse_reaction_name(line: str) -> Optional[str]:
    """
    Parse a reaction name from a subsection line.
    
    Examples:
        "Reaction 2.1.1 $   e + H(1s) \\rightarrow e + H(2p)$" -> "2.1.1"
        "Reaction 2.2.14 $  e + H_2^+(v) \\rightarrow H(1s) + H(n)$" -> "2.2.14"
    """
    match = re.search(r'Reaction\s+(\d+\.\d+\.\d+[ab]?)', line)
    if match:
        return match.group(1)
    return None


def _parse_coefficients(lines: list, start_idx: int) -> Tuple[np.ndarray, Dict, int]:
    """
    Parse polynomial coefficients from verbatim block.
    
    Returns:
        coeffs: Array of b0-b8 coefficients
        metadata: Dict with Tmin, sv_min, sv_max, error
        end_idx: Index after the parsed block
    """
    coeffs = np.zeros(9)
    metadata = {}
    
    idx = start_idx
    while idx < len(lines):
        line = lines[idx]
        
        # Look for coefficient lines (b0, b1, etc.)
        b_matches = re.findall(r'b(\d)\s+([+-]?\d+\.\d+e[+-]\d+)', line)
        for b_idx, b_val in b_matches:
            coeffs[int(b_idx)] = float(b_val)
        
        # Look for metadata line (Tmin, sv values)
        if 'Tmin' in line:
            tmin_match = re.search(r'Tmin\s+([+-]?\d+\.\d+e[+-]\d+)', line)
            if tmin_match:
                metadata['Tmin'] = float(tmin_match.group(1))
            
            sv_min_match = re.search(r'<sv>\(Tmin\)\s+([+-]?\d+\.\d+e[+-]\d+)', line)
            if sv_min_match:
                metadata['sv_min'] = float(sv_min_match.group(1))
            
            sv_max_match = re.search(r'<sv>max\s+([+-]?\d+\.\d+e[+-]\d+)', line)
            if sv_max_match:
                metadata['sv_max'] = float(sv_max_match.group(1))
            
            error_match = re.search(r'Error\s+([+-]?\d+\.\d+e[+-]\d+)', line)
            if error_match:
                metadata['error'] = float(error_match.group(1))
        
        # End of verbatim block
        if '\\end{verbatim}' in line:
            return coeffs, metadata, idx + 1
        
        idx += 1
    
    return coeffs, metadata, idx


def load_hydhel(filepath: Optional[str] = None) -> None:
    """
    Load and parse the HYDHEL database file.
    
    Parameters
    ----------
    filepath : str, optional
        Path to hydhel.tex file. If None, uses default location.
    """
    global _HYDHEL_COEFFICIENTS, _HYDHEL_METADATA, _HYDHEL_LOADED
    
    if filepath is None:
        filepath = _DEFAULT_HYDHEL_PATH
    
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"HYDHEL file not found: {filepath}")
    
    with open(filepath, 'r') as f:
        content = f.read()
    
    lines = content.split('\n')
    
    current_reaction = None
    in_subsection = False
    idx = 0
    
    while idx < len(lines):
        line = lines[idx]
        
        # Track when we enter a subsection
        if '\\subsection' in line:
            in_subsection = True
        
        # Look for reaction name inside or after subsection
        if in_subsection and 'Reaction' in line:
            reaction_name = _parse_reaction_name(line)
            if reaction_name:
                current_reaction = reaction_name
                in_subsection = False  # Reset after finding reaction
        
        # Look for start of verbatim block with coefficients
        if '\\begin{verbatim}' in line and current_reaction:
            coeffs, metadata, new_idx = _parse_coefficients(lines, idx + 1)
            # Only store if we actually found coefficients (b0 should be non-zero for valid reactions)
            if coeffs[0] != 0:
                _HYDHEL_COEFFICIENTS[current_reaction] = coeffs
                _HYDHEL_METADATA[current_reaction] = metadata
            idx = new_idx
            current_reaction = None  # Reset after processing
            continue
        
        idx += 1
    
    _HYDHEL_LOADED = True


def get_coefficients(reaction: str) -> np.ndarray:
    """
    Get polynomial coefficients for a reaction.
    
    Parameters
    ----------
    reaction : str
        Reaction identifier (e.g., "2.2.14", "2.1.5")
        
    Returns
    -------
    coeffs : np.ndarray
        Array of 9 polynomial coefficients (b0-b8)
    """
    if not _HYDHEL_LOADED:
        load_hydhel()
    
    if reaction not in _HYDHEL_COEFFICIENTS:
        raise KeyError(f"Reaction {reaction} not found in HYDHEL database. "
                      f"Available: {sorted(_HYDHEL_COEFFICIENTS.keys())}")
    
    return _HYDHEL_COEFFICIENTS[reaction].copy()


def get_metadata(reaction: str) -> Dict:
    """
    Get metadata for a reaction (Tmin, sv values, error).
    
    Parameters
    ----------
    reaction : str
        Reaction identifier (e.g., "2.2.14")
        
    Returns
    -------
    metadata : dict
        Dictionary with Tmin, sv_min, sv_max, error
    """
    if not _HYDHEL_LOADED:
        load_hydhel()
    
    if reaction not in _HYDHEL_METADATA:
        raise KeyError(f"Reaction {reaction} not found in HYDHEL database")
    
    return _HYDHEL_METADATA[reaction].copy()


# C++ Tmin thresholds for each reaction (returns k=0 below these temperatures)
# Reactions not listed here use global Tmin=0.1 clamp without returning 0
_REACTION_TMIN_THRESHOLDS = {
    "2.1.1": 0.4,
    "2.1.2": 0.4,
    "2.1.4a": 0.5,
    "2.1.5": 0.5,
    "2.2.2": 0.4,
    "2.2.3": 0.5,
    "2.2.4": 0.6,
    "2.2.5": 0.4,
    "2.2.6": 0.5,
    "2.2.7": 1.0,
    "2.2.8": 0.7,
    "2.2.9": 0.6,
    "2.2.10": 0.6,
    "2.2.11": 0.6,
    "2.2.13": 0.6,
    "2.2.16": 0.5,
}


def compute_rate(reaction: str, Te: np.ndarray) -> np.ndarray:
    """
    Compute rate coefficient using HYDHEL polynomial fit.
    
    k(T) = exp(sum_i b_i * (ln T)^i)
    
    This function matches the C++ behavior exactly:
    - Global clamp: T = max(min(T, 2e4), 0.1) (0.1 to 20000 eV range)
    - Per-reaction threshold: Returns k=0 if T < Tmin for certain reactions
    - Final check: k = max(k, 0)
    
    Parameters
    ----------
    reaction : str
        Reaction identifier (e.g., "2.2.14")
    Te : np.ndarray
        Temperature [eV]
        
    Returns
    -------
    k : np.ndarray
        Rate coefficient [cm³/s]
    """
    coeffs = get_coefficients(reaction)
    
    Te = np.atleast_1d(Te).astype(float)
    
    # C++ behavior: Global clamp to 0.1-20000 eV range (from reactionrates.cpp)
    Te_safe = np.clip(Te, 0.1, 2e4)
    
    # Compute rate using polynomial fit
    ln_Te = np.log(Te_safe)
    ln_k = np.zeros_like(Te_safe)
    
    for i, bi in enumerate(coeffs):
        ln_k += bi * np.power(ln_Te, i)
    
    k = np.exp(ln_k)
    
    # C++ behavior: Certain reactions return k=0 below their specific threshold
    if reaction in _REACTION_TMIN_THRESHOLDS:
        Tmin_threshold = _REACTION_TMIN_THRESHOLDS[reaction]
        k = np.where(Te < Tmin_threshold, 0.0, k)
    
    # C++ behavior: Ensure k >= 0
    k = np.maximum(k, 0.0)
    
    return k


def list_reactions() -> list:
    """List all available reactions in the database."""
    if not _HYDHEL_LOADED:
        load_hydhel()
    return sorted(_HYDHEL_COEFFICIENTS.keys())


# Auto-load on import if file exists
if _DEFAULT_HYDHEL_PATH.exists():
    try:
        load_hydhel()
    except Exception as e:
        print(f"Warning: Could not auto-load HYDHEL database: {e}")
