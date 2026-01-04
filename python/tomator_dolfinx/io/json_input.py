"""
JSON input file parser for Tomator simulations.

Parses the JSON input format used by C++ Tomator1D and extracts
parameters into a flat dictionary for the Python solver.

Supports two JSON formats:
1. Simple format: "param": value
2. Value-unit format: "param": {"value": value, "unit": "unit_string"}

The value-unit format is automatically detected and handled transparently.
"""

import json
from typing import Dict, Any, Union
from pathlib import Path


# Physical constants
KB = 1.380649e-23  # Boltzmann constant [J/K]


def get_value(data: Union[Dict, float, int, bool, str], default: Any = None) -> Any:
    """
    Extract value from a parameter that may be either:
    - A simple value (float, int, bool, str)
    - A dict with {"value": ..., "unit": ...} format
    
    Parameters
    ----------
    data : dict or scalar
        The parameter value, either plain or as value-unit dict.
    default : any, optional
        Default value if data is None.
        
    Returns
    -------
    value : scalar
        The extracted numeric/boolean/string value.
    """
    if data is None:
        return default
    if isinstance(data, dict):
        if 'value' in data:
            return data['value']
        # Empty dict or dict without 'value' key
        return default
    return data


def dict_get(d: Dict, key: str, default: Any = None) -> Any:
    """
    Get a value from a dict, handling value-unit format.
    
    Parameters
    ----------
    d : dict
        Dictionary to get value from.
    key : str
        Key to look up.
    default : any, optional
        Default value if key not found.
        
    Returns
    -------
    value : scalar
        The extracted value.
    """
    if key not in d:
        return default
    return get_value(d[key], default)


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

    params['operator_splitting'] = False  # Default to operator splitting
    
    # Magnetic field
    if 'magnetic_field' in raw:
        mf = raw['magnetic_field']
        params['Bt'] = dict_get(mf, 'Bt', 1.0)  # Toroidal field [T]
        params['Bv'] = dict_get(mf, 'Bv', 0.01)  # Vertical field [T]
        params['Bh'] = dict_get(mf, 'Bh', 0.001)  # Horizontal field [T]
    
    # Geometry - NOTE: C++ uses cm, we convert to m
    if 'toroidal_machine_geometry' in raw:
        geo = raw['toroidal_machine_geometry']
        params['R'] = dict_get(geo, 'R', 88.0) / 100.0  # Major radius [m]
        params['a'] = dict_get(geo, 'a', 28.0) / 100.0  # Minor radius [m]
        params['b'] = dict_get(geo, 'b', 75.0) / 100.0  # Vertical extent [m]
        params['lHFS'] = dict_get(geo, 'lHFS', 28.0) / 100.0  # HFS limiter distance [m]
        params['lLFS'] = dict_get(geo, 'lLFS', 28.0) / 100.0  # LFS limiter distance [m]
        params['nlimiters'] = dict_get(geo, 'nlimiters', 1)
        params['Vpl'] = dict_get(geo, 'Vpl', 1.0) / 1e6  # Plasma volume [m³]
    
    # Neutral pressure - convert from mbar to Pa (1 mbar = 100 Pa)
    if 'neutral_pressure' in raw:
        np_data = raw['neutral_pressure']
        params['pHe'] = dict_get(np_data, 'pHe', 0.0) * 100.0  # He pressure [Pa]
        params['pH2'] = dict_get(np_data, 'pH2', 0.0) * 100.0  # H2 pressure [Pa]
    
    # RF power
    if 'rf_power' in raw:
        rf = raw['rf_power']
        params['Prf'] = dict_get(rf, 'Prf', 0.0)  # RF power [kW]
        params['freq'] = dict_get(rf, 'freq', 1e6)  # Frequency [Hz]
        params['dtpramp'] = dict_get(rf, 'dtpramp', 0.0)  # Power ramp time [s]
    
    # Coupled power mode flags
    if 'type' in raw:
        ptype = raw['type']
        params['bnefix'] = dict_get(ptype, 'bnefix', False)
        params['bfixpowerfrac'] = dict_get(ptype, 'bfixpowerfrac', False)
        params['bproptone'] = dict_get(ptype, 'bproptone', False)
        params['bnopower'] = dict_get(ptype, 'bnopower', False)
    
    # General EC parameters
    if 'general_ec' in raw:
        ec = raw['general_ec']
        params['Rdep'] = dict_get(ec, 'Rdep', 88.0) / 100.0  # Deposition radius [m]
        params['pecabs0'] = dict_get(ec, 'pecabs0', 0.1)      # Initial absorbed fraction
        params['widthech'] = dict_get(ec, 'widthech', 2.0) / 100.0  # EC width [m]
        params['echbackground'] = dict_get(ec, 'echbackground', 1e-5)  # Background fraction
    
    # necfix mode parameters
    if 'necfix' in raw:
        nf = raw['necfix']
        params['ic'] = dict_get(nf, 'ic', 90)                # Control index
        params['necfix'] = dict_get(nf, 'necfix', 1e13) * 1e6  # Target density [m^-3]
        params['P_KP'] = dict_get(nf, 'P_KP', 10.0)          # PID proportional gain
        params['P_KI'] = dict_get(nf, 'P_KI', 10000.0)       # PID integral gain
        params['P_KD'] = dict_get(nf, 'P_KD', 0.1)           # PID derivative gain
        params['P_KI_ini'] = dict_get(nf, 'P_KI_ini', 0.0)   # Initial integral term value
    
    # Physics flags
    if 'physics_to_include' in raw:
        phys = raw['physics_to_include']
        params['bH'] = dict_get(phys, 'bH', True)
        params['bH2'] = dict_get(phys, 'bH2', False)
        params['bHe'] = dict_get(phys, 'bHe', False)
        params['bion'] = dict_get(phys, 'bion', True)
        params['bcx'] = dict_get(phys, 'bcx', True)
        params['belas'] = dict_get(phys, 'belas', False)
        params['bcoulomb'] = dict_get(phys, 'bcoulomb', False)
        params['bimpur'] = dict_get(phys, 'bimpur', False)
        params['btranspions'] = dict_get(phys, 'btranspions', True)
        params['btranspneut'] = dict_get(phys, 'btranspneut', True)
    
    # Diffusion parameters
    if 'diffusion' in raw:
        diff = raw['diffusion']
        
        # Support new nested structure: {"Dfix": {"bool": true, "Dfix": 3330.0, "unit": "cm²/s"}}
        # as well as old flat structure: {"bDfix": {"value": true}, "Dfix": {"value": 3330.0}}
        
        # Check for new nested structure first
        if 'Dfix' in diff and isinstance(diff['Dfix'], dict) and 'bool' in diff['Dfix']:
            # New structure
            params['bDfix'] = diff['Dfix'].get('bool', False)
            params['Dfix'] = diff['Dfix'].get('Dfix', 1.0) / 1e4  # Convert cm²/s to m²/s
        else:
            # Old structure
            params['bDfix'] = dict_get(diff, 'bDfix', True)
            params['Dfix'] = dict_get(diff, 'Dfix', 1.0) / 1e4  # Convert cm²/s to m²/s
        
        if 'Dbohm' in diff and isinstance(diff['Dbohm'], dict) and 'bool' in diff['Dbohm']:
            # New structure
            params['bDbohm'] = diff['Dbohm'].get('bool', False)
            params['Dfact'] = diff['Dbohm'].get('Dfact', 1.0)
        else:
            # Old structure
            params['bDbohm'] = dict_get(diff, 'bDbohm', False)
            params['Dfact'] = dict_get(diff, 'Dfact', 1.0)
        
        # Support both Dgyrogeom (new name) and Dscaling (legacy name)
        if 'Dgyrogeom' in diff and isinstance(diff['Dgyrogeom'], dict) and 'bool' in diff['Dgyrogeom']:
            # New structure with Dgyrogeom
            params['bDgyrogeom'] = diff['Dgyrogeom'].get('bool', False)
            if params['bDgyrogeom']:
                params['Dfact'] = diff['Dgyrogeom'].get('Dfact', params.get('Dfact', 1.0))
        elif 'Dscaling' in diff and isinstance(diff['Dscaling'], dict) and 'bool' in diff['Dscaling']:
            # Legacy structure with Dscaling
            params['bDgyrogeom'] = diff['Dscaling'].get('bool', False)
            if params['bDgyrogeom']:
                params['Dfact'] = diff['Dscaling'].get('Dfact', params.get('Dfact', 1.0))
        else:
            # Old flat structure
            params['bDgyrogeom'] = dict_get(diff, 'bDscaling', dict_get(diff, 'bDgyrogeom', False))
        
        # Set diffusion type based on priority: gyrogeom > bohm > fixed
        if params.get('bDgyrogeom', False):
            params['diffusion_type'] = 'gyrogeom'
        elif params.get('bDbohm', False):
            params['diffusion_type'] = 'bohm'
        elif params.get('bDfix', False):
            params['diffusion_type'] = 'fixed'
        else:
            params['diffusion_type'] = 'fixed'
    
    # Advection parameters (also accepts old 'convection' key for backward compatibility)
    adv_key = 'advection' if 'advection' in raw else 'convection'
    if adv_key in raw:
        adv = raw[adv_key]
        
        # Support new nested structure: {"Vfix": {"bool": true, "value": 333.0, "unit": "cm/s"}}
        # as well as old flat structure: {"bVfix": {"value": true}, "Vfix": {"value": 333.0}}
        
        # Check for new nested structure first
        if 'Vfix' in adv and isinstance(adv['Vfix'], dict) and 'bool' in adv['Vfix']:
            # New structure
            params['bVfix'] = adv['Vfix'].get('bool', False)
            params['Vfix'] = adv['Vfix'].get('value', 0.0) / 100.0  # Convert cm/s to m/s
        else:
            # Old structure
            params['bVfix'] = dict_get(adv, 'bVfix', True)
            params['Vfix'] = dict_get(adv, 'Vfix', 0.0) / 100.0  # Convert cm/s to m/s
        
        if 'Vscaling' in adv and isinstance(adv['Vscaling'], dict) and 'bool' in adv['Vscaling']:
            # New structure
            params['bVscaling'] = adv['Vscaling'].get('bool', False)
            params['veq'] = adv['Vscaling'].get('veq', 8)
            params['Vfact'] = adv['Vscaling'].get('Vfact', 1.0)
        else:
            # Old structure
            params['bVscaling'] = dict_get(adv, 'bVscaling', False)
            params['veq'] = dict_get(adv, 'veq', 8)
            params['Vfact'] = dict_get(adv, 'Vfact', 1.0)
        
        # Set advection type
        if params['bVfix']:
            params['advection_type'] = 'fixed'
        else:
            params['advection_type'] = 'pressure'
    
    # Transport dict for TransportManager
    # Include all transport-related flags for proper model selection
    params['transport'] = {
        'Dfix': params.get('Dfix', 0.01),
        'Vfix': params.get('Vfix', 0.0),
        'diffusion_type': params.get('diffusion_type', 'fixed'),
        'advection_type': params.get('advection_type', 'fixed'),
        # Include bool flags for model selection in initialize_from_params
        'bDfix': params.get('bDfix', True),
        'bDbohm': params.get('bDbohm', False),
        'bDgyrogeom': params.get('bDgyrogeom', False),
        'Dfact': params.get('Dfact', 1.0),  # Scaling factor (Dfsave) for Bohm/Gyrogeom
        'bVfix': params.get('bVfix', True),
        'bVscaling': params.get('bVscaling', False),
        'veq': params.get('veq', 8),
        'Vfact': params.get('Vfact', 1.0),
    }
    
    # Initial conditions
    if 'initial_conditions' in raw:
        ic = raw['initial_conditions']
        # Temperatures [eV]
        params['Te0'] = dict_get(ic, 'Te0', 10.0)
        params['Ta0'] = dict_get(ic, 'Ta0', 0.026)  # Room temperature / neutral temp
        
        # Vacuum density floor (critical for stability)
        # C++ uses cm^-3, we keep it in cm^-3 and convert in solver
        params['nevac'] = dict_get(ic, 'nevac', 1.0)  # [cm^-3]
        
        # Densities - convert from cm^-3 to m^-3
        params['nH0'] = dict_get(ic, 'nH0', 1e10) * 1e6
        params['nHi0'] = dict_get(ic, 'nHi0', 1e10) * 1e6
        params['nH20'] = dict_get(ic, 'nH20', 1e10) * 1e6
        params['nH2i0'] = dict_get(ic, 'nH2i0', 1e8) * 1e6
        params['nH3i0'] = dict_get(ic, 'nH3i0', 1e8) * 1e6
        params['nHeII0'] = dict_get(ic, 'nHeII0', 1e10) * 1e6
        params['nHeIII0'] = dict_get(ic, 'nHeIII0', 1e6) * 1e6
        
        # nHeI0 and nH20: compute from pressure if not explicitly given
        # This matches C++ Tomator1D behavior: nHeI0 = computeN(pHe, Ta0)
        Ta0 = params.get('Ta0', 0.026)  # Neutral temperature [eV]
        
        if 'nHeI0' in ic:
            params['nHeI0'] = dict_get(ic, 'nHeI0', 1e10) * 1e6  # Convert cm^-3 to m^-3
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
        params['rmaxini'] = dict_get(ic, 'rmaxini', params.get('R', 0.88) * 100) / 100.0  # [m]
        params['widthini'] = dict_get(ic, 'widthini', 5.0) / 100.0  # [m]
        params['nebackgroundl'] = dict_get(ic, 'nebackgroundl', 1e-5)  # HFS background fraction
        params['nebackgroundr'] = dict_get(ic, 'nebackgroundr', 1e-3)  # LFS background fraction
        
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
        params['RH'] = dict_get(edge, 'RH', 0.5)  # Reflection coefficient for H
        params['REH'] = dict_get(edge, 'REH', 0.9)  # Energy reflection coefficient for neutrals
        params['gEd'] = dict_get(edge, 'gEd', 5/3)  # Energy flux factor for diffusion
        params['gEv'] = dict_get(edge, 'gEv', 5/3)  # Energy flux factor for advection
        params['gEdn'] = dict_get(edge, 'gEdn', 5/3)  # Energy flux factor for neutral diffusion
        params['gEe'] = dict_get(edge, 'gEe', 5/3)  # Energy flux factor for electrons
        
        # Fixed BC option for ions (C++ fixBCs)
        # If true, use fixed decay length instead of physics-based calculation
        if 'fixBC' in edge and isinstance(edge['fixBC'], dict) and 'bool' in edge['fixBC']:
            params['fixBC'] = edge['fixBC'].get('bool', False)
            params['fixBC_value'] = edge['fixBC'].get('value', 2.0) / 100.0  # Convert cm to m
        else:
            params['fixBC'] = dict_get(edge, 'fixBC', False)
            params['fixBC_value'] = 0.02  # Default 2 cm in m
    
    # Decay lengths - initial values only, updated dynamically based on actual D
    params['decay_length_hfs'] = 0.01  # [m] - will be recomputed from physics
    params['decay_length_lfs'] = 0.01  # [m] - will be recomputed from physics
    
    # Simulation grid
    if 'simulation_grid' in raw:
        grid = raw['simulation_grid']
        params['nmeshp'] = dict_get(grid, 'nmeshp', 101)
        # FEM polynomial degree: 1=linear, 2=quadratic, 3=cubic
        # Higher degree gives wider stencil and can improve stability
        params['fem_degree'] = dict_get(grid, 'fem_degree', 1)
    
    # Time stepping
    if 'time_step' in raw:
        ts = raw['time_step']
        params['t0'] = dict_get(ts, 't0', 0.0)
        params['tmainend'] = dict_get(ts, 'tmainend', 1e-3)
        params['accur'] = dict_get(ts, 'accur', 0.05)
        params['dtmax'] = dict_get(ts, 'dtmax', 1e-5)
        params['dtmin'] = dict_get(ts, 'dtmin', 1e-10)
        params['dtinit'] = dict_get(ts, 'dtinit', 1e-9)
    
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
        params['Nlog'] = dict_get(out, 'Nlog', 100)
        params['dtsave'] = dict_get(out, 'dtsave', 1e-4)
        params['output_interval'] = dict_get(out, 'dtsave', 1e-4)
        # Profiling flag (can be plain bool or {"value": bool})
        if 'profile' in out:
            profile_val = out['profile']
            params['profile'] = profile_val if isinstance(profile_val, bool) else profile_val.get('value', False)
    
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
            'advection_type': 'fixed',
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
