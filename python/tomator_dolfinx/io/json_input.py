"""
JSON input file parser for Tomator simulations.

Parses the JSON input format used by C++ Tomator1D and extracts
parameters into a flat dictionary for the Python solver.
"""

import json
from typing import Dict, Any
from pathlib import Path


# Physical constants
KB = 1.380649e-23  # Boltzmann constant [J/K]


def compute_density_from_pressure(p_pa: float, T_eV: float) -> float:
    """
    Compute number density from pressure using ideal gas law.
    n = p / (kB * T)
    
    Parameters
    ----------
    p_pa : float
        Pressure in Pascal [Pa].
    T_eV : float
        Temperature in electronvolts [eV].
        
    Returns
    -------
    n : float
        Number density in m^-3.
    """
    T_K = T_eV * 11600.0  # Convert eV to K
    return p_pa / (KB * T_K)  # [m^-3]


def load_input_file(filename: str) -> Dict[str, Any]:
    """
    Load and parse a Tomator JSON input file.
    
    Parameters
    ----------
    filename : str
        Path to JSON input file.
        
    Returns
    -------
    params : dict
        Flattened parameter dictionary.
    """
    with open(filename, 'r') as f:
        raw = json.load(f)
    
    # Flatten nested structure into single dict
    params = {}
    
    # Magnetic field
    if 'magnetic_field' in raw:
        mf = raw['magnetic_field']
        params['Bt'] = mf.get('Bt', 1.0)  # Toroidal field [T]
        params['Bv'] = mf.get('Bv', 0.01)  # Vertical field [T]
        params['Bh'] = mf.get('Bh', 0.001)  # Horizontal field [T]
    
    # Geometry - NOTE: C++ uses cm, we convert to m
    if 'toroidal_machine_geometry' in raw:
        geo = raw['toroidal_machine_geometry']
        params['R'] = geo.get('R', 88.0) / 100.0  # Major radius [m]
        params['a'] = geo.get('a', 28.0) / 100.0  # Minor radius [m]
        params['b'] = geo.get('b', 75.0) / 100.0  # Vertical extent [m]
        params['lHFS'] = geo.get('lHFS', 28.0) / 100.0  # HFS limiter distance [m]
        params['lLFS'] = geo.get('lLFS', 28.0) / 100.0  # LFS limiter distance [m]
        params['nlimiters'] = geo.get('nlimiters', 1)
        params['Vpl'] = geo.get('Vpl', 1.0) / 1e6  # Plasma volume [m³]
    
    # Neutral pressure - convert from mbar to Pa (1 mbar = 100 Pa)
    if 'neutral_pressure' in raw:
        np_data = raw['neutral_pressure']
        params['pHe'] = np_data.get('pHe', 0.0) * 100.0  # He pressure [Pa]
        params['pH2'] = np_data.get('pH2', 0.0) * 100.0  # H2 pressure [Pa]
    
    # RF power
    if 'rf_power' in raw:
        rf = raw['rf_power']
        params['Prf'] = rf.get('Prf', 0.0)  # RF power [kW]
        params['freq'] = rf.get('freq', 1e6)  # Frequency [Hz]
    
    # Physics flags
    if 'physics_to_include' in raw:
        phys = raw['physics_to_include']
        params['bH'] = phys.get('bH', True)
        params['bH2'] = phys.get('bH2', False)
        params['bHe'] = phys.get('bHe', False)
        params['bion'] = phys.get('bion', True)
        params['bcx'] = phys.get('bcx', True)
        params['belas'] = phys.get('belas', False)
        params['bcoulomb'] = phys.get('bcoulomb', False)
        params['bimpur'] = phys.get('bimpur', False)
        params['btranspions'] = phys.get('btranspions', True)
        params['btranspneut'] = phys.get('btranspneut', True)
    
    # Diffusion parameters
    if 'diffusion' in raw:
        diff = raw['diffusion']
        params['bDfix'] = diff.get('bDfix', True)
        params['Dfix'] = diff.get('Dfix', 1.0) / 1e4  # Convert cm²/s to m²/s
        params['bDbohm'] = diff.get('bDbohm', False)
        params['Dfact'] = diff.get('Dfact', 1.0)
        
        # Set diffusion type
        if diff.get('bDfix', True):
            params['diffusion_type'] = 'fixed'
        elif diff.get('bDbohm', False):
            params['diffusion_type'] = 'bohm'
        else:
            params['diffusion_type'] = 'classical'
    
    # Convection parameters
    if 'convection' in raw:
        conv = raw['convection']
        params['bVfix'] = conv.get('bVfix', True)
        params['Vfix'] = conv.get('Vfix', 0.0) / 100.0  # Convert cm/s to m/s
        params['Vfact'] = conv.get('Vfact', 1.0)
        
        # Set convection type
        if conv.get('bVfix', True):
            params['convection_type'] = 'fixed'
        else:
            params['convection_type'] = 'pressure'
    
    # Transport dict for TransportManager
    params['transport'] = {
        'Dfix': params.get('Dfix', 0.01),
        'Vfix': params.get('Vfix', 0.0),
        'diffusion_type': params.get('diffusion_type', 'fixed'),
        'convection_type': params.get('convection_type', 'fixed'),
    }
    
    # Initial conditions
    if 'initial_conditions' in raw:
        ic = raw['initial_conditions']
        # Temperatures [eV]
        params['Te0'] = ic.get('Te0', 10.0)
        params['Ta0'] = ic.get('Ta0', 0.026)  # Room temperature / neutral temp
        
        # Vacuum density floor (critical for stability)
        # C++ uses cm^-3, we keep it in cm^-3 and convert in solver
        params['nevac'] = ic.get('nevac', 1.0)  # [cm^-3]
        
        # Densities - convert from cm^-3 to m^-3
        params['nH0'] = ic.get('nH0', 1e10) * 1e6
        params['nHi0'] = ic.get('nHi0', 1e10) * 1e6
        params['nH20'] = ic.get('nH20', 1e10) * 1e6
        params['nH2i0'] = ic.get('nH2i0', 1e8) * 1e6
        params['nH3i0'] = ic.get('nH3i0', 1e8) * 1e6
        params['nHeII0'] = ic.get('nHeII0', 1e10) * 1e6
        params['nHeIII0'] = ic.get('nHeIII0', 1e6) * 1e6
        
        # nHeI0 and nH20: compute from pressure if not explicitly given
        # This matches C++ Tomator1D behavior: nHeI0 = computeN(pHe, Ta0)
        Ta0 = params.get('Ta0', 0.026)  # Neutral temperature [eV]
        
        if 'nHeI0' in ic:
            params['nHeI0'] = ic['nHeI0'] * 1e6  # Convert cm^-3 to m^-3
        elif 'pHe' in params and params['pHe'] > 0:
            # Compute from He pressure [Pa] and neutral temp [eV]
            params['nHeI0'] = compute_density_from_pressure(params['pHe'], Ta0)
        else:
            params['nHeI0'] = 1e16  # Default [m^-3]
        
        # Set nHeI boundary condition = initial value (fixed neutral density at edges)
        params['nHeI_bc'] = params['nHeI0']
            
        if 'nH20' not in ic and 'pH2' in params and params['pH2'] > 0:
            # Compute from H2 pressure [Pa] and neutral temp [eV]
            params['nH20'] = compute_density_from_pressure(params['pH2'], Ta0)
        
        # Set nH2 boundary condition = initial value
        params['nH2_bc'] = params.get('nH20', 1e15)
        
        # Profile parameters - for Gaussian initialization
        params['rmaxini'] = ic.get('rmaxini', params.get('R', 0.88) * 100) / 100.0  # [m]
        params['widthini'] = ic.get('widthini', 5.0) / 100.0  # [m]
        params['nebackgroundl'] = ic.get('nebackgroundl', 1e-5)  # HFS background fraction
        params['nebackgroundr'] = ic.get('nebackgroundr', 1e-3)  # LFS background fraction
        
    params['initial_conditions'] = {
        'Te0': params.get('Te0', 10.0),
        'Ta0': params.get('Ta0', 0.026),
        'TH0': params.get('Ta0', 0.026),
        'THi0': params.get('Te0', 10.0),
        'nH0': params.get('nH0', 1e16),
        'nHi0': params.get('nHi0', 1e16),
        'nH20': params.get('nH20', 1e15),
        'nH2i0': params.get('nH2i0', 1e12),
        'nH3i0': params.get('nH3i0', 1e12),
        'nHeI0': params.get('nHeI0', 1e15),
        'nHeII0': params.get('nHeII0', 1e14),
        'nHeIII0': params.get('nHeIII0', 1e12),
        'nevac': params.get('nevac', 1.0),  # cm^-3
        # Gaussian profile parameters
        'rmaxini': params.get('rmaxini', 0.92),  # [m]
        'widthini': params.get('widthini', 0.05),  # [m]
        'nebackgroundl': params.get('nebackgroundl', 1e-5),
        'nebackgroundr': params.get('nebackgroundr', 1e-3),
    }
    
    # Edge/boundary conditions
    if 'edge_conditions' in raw:
        edge = raw['edge_conditions']
        params['RH'] = edge.get('RH', 0.5)  # Reflection coefficient for H
        params['gEd'] = edge.get('gEd', 5/3)  # Energy flux factor for diffusion
        params['gEv'] = edge.get('gEv', 5/3)  # Energy flux factor for convection
    
    # Decay lengths (derive from geometry if not specified)
    params['decay_length_hfs'] = params.get('lHFS', 0.2) * 0.1  # ~10% of limiter distance
    params['decay_length_lfs'] = params.get('lLFS', 0.2) * 0.1
    
    # Simulation grid
    if 'simulation_grid' in raw:
        grid = raw['simulation_grid']
        params['nmeshp'] = grid.get('nmeshp', 101)
    
    # Time stepping
    if 'time_step' in raw:
        ts = raw['time_step']
        params['t0'] = ts.get('t0', 0.0)
        params['tmainend'] = ts.get('tmainend', 1e-3)
        params['accur'] = ts.get('accur', 0.05)
        params['dtmax'] = ts.get('dtmax', 1e-5)
        params['dtmin'] = ts.get('dtmin', 1e-10)
        params['dtinit'] = ts.get('dtinit', 1e-9)
    
    params['time_step'] = {
        'dtinit': params.get('dtinit', 1e-9),
        'dtmin': params.get('dtmin', 1e-10),
        'dtmax': params.get('dtmax', 1e-5),
        'accur': params.get('accur', 0.05),
        'decay_length_hfs': params.get('decay_length_hfs', 0.02),
        'decay_length_lfs': params.get('decay_length_lfs', 0.02),
        'nevac': params.get('nevac', 1.0),  # Vacuum density floor [cm^-3]
    }
    
    # Output parameters
    if 'output_parameters' in raw:
        out = raw['output_parameters']
        params['Nlog'] = out.get('Nlog', 100)
        params['dtsave'] = out.get('dtsave', 1e-4)
        params['output_interval'] = out.get('dtsave', 1e-4)
    
    return params


def load_grid_file(filename: str) -> list:
    """
    Load radial grid from JSON file.
    
    Expected format:
    {
        "radial_positions": [r0, r1, r2, ..., rN]
    }
    
    or simply an array of positions.
    
    Parameters
    ----------
    filename : str
        Path to grid JSON file.
        
    Returns
    -------
    positions : list
        Radial positions [m].
    """
    with open(filename, 'r') as f:
        data = json.load(f)
    
    if isinstance(data, list):
        return data
    elif 'radial_positions' in data:
        return data['radial_positions']
    else:
        raise ValueError(f"Invalid grid file format: {filename}")


def create_default_params() -> Dict[str, Any]:
    """
    Create default parameter dictionary.
    
    Returns
    -------
    params : dict
        Default parameters.
    """
    return {
        # Geometry
        'R': 0.88,
        'a': 0.28,
        'lHFS': 0.28,
        'lLFS': 0.28,
        
        # Magnetic field
        'Bt': 1.5,
        'Bv': 0.01,
        'Bh': 0.001,
        
        # Physics
        'bH': True,
        'bH2': False,
        'bHe': False,
        'bcx': True,
        
        # Transport
        'transport': {
            'Dfix': 0.1,
            'Vfix': 0.0,
            'diffusion_type': 'fixed',
            'convection_type': 'fixed',
        },
        
        # Initial conditions
        'initial_conditions': {
            'Te0': 10.0,
            'TH0': 0.026,
            'THi0': 10.0,
            'nH0': 1e16,
            'nHi0': 1e18,
        },
        
        # Time stepping
        'time_step': {
            'dtinit': 1e-9,
            'dtmin': 1e-10,
            'dtmax': 1e-5,
            'accur': 0.05,
            'decay_length_hfs': 0.02,
            'decay_length_lfs': 0.02,
        },
        
        # Grid
        'nmeshp': 101,
        
        # Output
        'tmainend': 1e-3,
        'output_interval': 1e-5,
    }
