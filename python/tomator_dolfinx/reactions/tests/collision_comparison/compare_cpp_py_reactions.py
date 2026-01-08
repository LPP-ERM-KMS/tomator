#!/usr/bin/env python3
"""
Compare dn and dE between C++ collisions.cpp and Python collisions.py
for each individual reaction flag.

This script:
1. Runs Python collisions with a single reaction flag enabled
2. Runs C++ test program with the same flag (compiled against actual Tomator source)
3. Compares the results and reports differences

Usage:
    python compare_cpp_py_reactions.py [--flag FLAG_NAME] [--all] [--tolerance TOL]
    
To compile the C++ test program manually:
    cd /home/ITER/wautert/Documents/Tomator-RMA/tomator/src/build
    g++ -O2 -std=c++17 -fopenmp -I.. -I../Eigen \
        ../tests/test_single_collision.cpp \
        ../Funcs/*.cpp ../Vars/*.cpp ../Vars/*.c \
        -o test_single_collision
"""

import numpy as np
import subprocess
import json
import os
import sys
import argparse

sys.path.insert(0, '/home/ITER/wautert/Documents/Tomator-RMA/tomator/python')

from tomator_dolfinx.reactions import collisions

# Paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TOMATOR_SRC = '/home/ITER/wautert/Documents/Tomator-RMA/tomator/src'
CPP_TEST_PROG = os.path.join(TOMATOR_SRC, 'build', 'test_single_collision')
CPP_SOURCE = os.path.join(TOMATOR_SRC, 'tests', 'test_single_collision.cpp')
BENCHMARK_PYTHON_DIR = '/home/ITER/wautert/Documents/Tomator-RMA/tomator/python/benchmark/python'

# Unit conversions
CM3_TO_M3 = 1e-6
M3_TO_CM3 = 1e6

# Test conditions (SI units for Python, C++ uses CGS internally)
# These values are in m^-3 for Python, C++ test program expects cm^-3
TEST_CONDITIONS = {
    'Te': 10.0,       # eV
    'ne': 1e18,       # m^-3 (= 1e12 cm^-3)
    'TH': 0.5,        # eV
    'nH': 1e16,       # m^-3 (= 1e10 cm^-3)
    'THi': 8.0,       # eV
    'nHi': 9e17,      # m^-3 (= 9e11 cm^-3)
    'TH2': 0.3,       # eV
    'nH2': 5e15,      # m^-3 (= 5e9 cm^-3)
    'TH2i': 2.0,      # eV
    'nH2i': 1e15,     # m^-3 (= 1e9 cm^-3)
    'TH3i': 1.5,      # eV
    'nH3i': 1e14,     # m^-3 (= 1e8 cm^-3)
    'THeI': 0.5,      # eV
    'nHeI': 1e16,     # m^-3 (= 1e10 cm^-3)
    'THeII': 5.0,     # eV
    'nHeII': 5e15,    # m^-3 (= 5e9 cm^-3)
    'THeIII': 10.0,   # eV
    'nHeIII': 1e14,   # m^-3 (= 1e8 cm^-3)
    # Disable rate limiting for fair comparison with C++ (which doesn't limit)
    'dt': None,
    'accur': None,
}

# All individual reaction flags
ALL_FLAGS = {
    'bH': ['bH_exc', 'bH_ion', 'bH_3body', 'bH_rec'],
    'bH2': ['bH2_elas', 'bH2_exc', 'bH2_dis', 'bH2_ion', 'bH2i_rec', 'bH2i_disrec', 
            'bH2_dision', 'bH3i_disrec', 'bH2i_dis', 'bH2i_disexc', 'bH3i_dis'],
    'bHe': ['bHeI_ion', 'bHeII_ion', 'bHe_cooling', 'bHeII_rec', 'bHeIII_rec'],
    'bcx': ['bcx_HiH', 'bcx_HiH2', 'bcx_H2iH2', 'bcx_HeIIH', 'bcx_HeIIHeI', 
            'bcx_HeIIIH', 'bcx_HeIIIHeI'],
    'bion': ['bion_HiH_exca', 'bion_HiH_excb', 'bion_HiH_ion', 'bion_HiH2_exca',
             'bion_HiH2_excb', 'bion_HiH2_325', 'bion_HiH2i_326', 'bion_H2iH2_H3i',
             'bion_HiHeI_ion', 'bion_HeIIH2_cxdis'],
    'belas': ['belas_HiH', 'belas_HH', 'belas_HH2', 'belas_HiH2', 'belas_H2iH',
              'belas_H2iH2', 'belas_H3iH', 'belas_H3iH2', 'belas_H2H2', 'belas_HHeI',
              'belas_HeIIH', 'belas_HiHeI', 'belas_HeIIHeI', 'belas_HeIHeI',
              'belas_HeIH2', 'belas_HeIIH2'],
    'bcoulomb': ['COULOMB_ALL'],  # Coulomb collisions tested as a group
}

def get_all_false_flags():
    """Return all individual flags set to False (excluding special flags)."""
    flags = {}
    for section_flags in ALL_FLAGS.values():
        for flag in section_flags:
            # Skip special flags that aren't actual Python kwargs
            if flag in ['COULOMB_ALL']:
                continue
            flags[flag] = False
    return flags


def get_all_true_flags():
    """Return all individual flags set to True (excluding special flags)."""
    flags = {}
    for section_flags in ALL_FLAGS.values():
        for flag in section_flags:
            # Skip special flags that aren't actual Python kwargs
            if flag in ['COULOMB_ALL']:
                continue
            flags[flag] = True
    return flags


def run_python_collision(flag_name, conditions=None):
    """Run Python collision with a single flag enabled (or all flags if flag_name='ALL_ENABLED')."""
    if conditions is None:
        conditions = TEST_CONDITIONS
    tc = conditions
    
    # Handle special cases
    if flag_name == 'ALL_ENABLED':
        flags = get_all_true_flags()
        include_coulomb = True  # Include Coulomb when running all
    elif flag_name == 'COULOMB_ALL':
        # Coulomb only - disable all other flags
        flags = get_all_false_flags()
        include_coulomb = True
    else:
        # All individual flags False except the one we're testing
        flags = get_all_false_flags()
        flags[flag_name] = True
        include_coulomb = False
    
    # Convert to numpy arrays
    result = collisions.compute_collision_sources(
        ne=np.array([tc['ne']]), Te=np.array([tc['Te']]),
        nH=np.array([tc['nH']]), TH=np.array([tc['TH']]),
        nHi=np.array([tc['nHi']]), THi=np.array([tc['THi']]),
        nH2=np.array([tc['nH2']]), TH2=np.array([tc['TH2']]),
        nH2i=np.array([tc['nH2i']]), TH2i=np.array([tc['TH2i']]),
        nH3i=np.array([tc['nH3i']]), TH3i=np.array([tc['TH3i']]),
        nHeI=np.array([tc['nHeI']]), THeI=np.array([tc['THeI']]),
        nHeII=np.array([tc['nHeII']]), THeII=np.array([tc['THeII']]),
        nHeIII=np.array([tc['nHeIII']]), THeIII=np.array([tc['THeIII']]),
        # Enable all section flags - individual flags control reactions
        include_H=True, include_H2=True, include_He=True,
        include_cx=True, include_ion=True, include_elastic=True,
        include_coulomb=include_coulomb,
        **flags
    )
    
    dn, dE, nu = result
    
    # Convert to dict with scalar values
    dn_dict = {k: float(v[0]) for k, v in dn.items()}
    dE_dict = {k: float(v[0]) for k, v in dE.items()}
    nu_dict = {k: float(v[0]) for k, v in nu.items()}
    
    return dn_dict, dE_dict, nu_dict


def run_cpp_collision(flag_name, conditions=None):
    """
    Run C++ collision test with a single flag enabled.
    Returns dn, dE dicts parsed from C++ output.
    C++ outputs in CGS (cm^-3/s), we convert to SI (m^-3/s).
    
    If conditions is provided, passes them to C++ via environment variables.
    Conditions should be in SI units (m^-3 for densities) - will be converted to CGS for C++.
    """
    if conditions is None:
        conditions = TEST_CONDITIONS
        
    if not os.path.exists(CPP_TEST_PROG):
        print(f"C++ test program not found: {CPP_TEST_PROG}")
        print("Attempting to compile...")
        compile_cpp()
        
    if not os.path.exists(CPP_TEST_PROG):
        print("Failed to compile C++ test program")
        return None, None, None
    
    # Build environment variables for C++ (convert SI to CGS for densities)
    env = os.environ.copy()
    env['TEST_Te'] = str(conditions['Te'])
    env['TEST_ne'] = str(conditions['ne'] * 1e-6)  # m^-3 to cm^-3
    env['TEST_TH'] = str(conditions['TH'])
    env['TEST_nH'] = str(conditions['nH'] * 1e-6)
    env['TEST_THi'] = str(conditions['THi'])
    env['TEST_nHi'] = str(conditions['nHi'] * 1e-6)
    env['TEST_TH2'] = str(conditions['TH2'])
    env['TEST_nH2'] = str(conditions['nH2'] * 1e-6)
    env['TEST_TH2i'] = str(conditions['TH2i'])
    env['TEST_nH2i'] = str(conditions['nH2i'] * 1e-6)
    env['TEST_TH3i'] = str(conditions['TH3i'])
    env['TEST_nH3i'] = str(conditions['nH3i'] * 1e-6)
    env['TEST_THeI'] = str(conditions['THeI'])
    env['TEST_nHeI'] = str(conditions['nHeI'] * 1e-6)
    env['TEST_THeII'] = str(conditions['THeII'])
    env['TEST_nHeII'] = str(conditions['nHeII'] * 1e-6)
    env['TEST_THeIII'] = str(conditions['THeIII'])
    env['TEST_nHeIII'] = str(conditions['nHeIII'] * 1e-6)
    
    # Run C++ program with flag name as argument
    try:
        result = subprocess.run(
            [CPP_TEST_PROG, flag_name],
            capture_output=True,
            text=True,
            timeout=10,
            env=env
        )
        
        if result.returncode != 0:
            print(f"C++ program error: {result.stderr}")
            return None, None, None
        
        # Parse JSON output from C++
        # The C++ program may print initialization messages before JSON
        # Find the JSON part starting with '{'
        stdout = result.stdout
        json_start = stdout.find('{')
        if json_start == -1:
            print(f"No JSON found in C++ output: {stdout}")
            return None, None, None
        json_str = stdout[json_start:]
        
        output = json.loads(json_str)
        
        # Convert CGS to SI: multiply by 1e6 (cm^-3 -> m^-3)
        dn_cgs = output.get('dn', {})
        dE_cgs = output.get('dE', {})
        nu_cgs = output.get('nu', {})
        
        # Convert dn from cm^-3/s to m^-3/s (multiply by 1e6)
        dn_si = {k: v * 1e6 for k, v in dn_cgs.items()}
        # Convert dE from eV*cm^-3/s to eV*m^-3/s (multiply by 1e6)
        dE_si = {k: v * 1e6 for k, v in dE_cgs.items()}
        # nu is in 1/s, no conversion needed
        nu_si = nu_cgs
        
        return dn_si, dE_si, nu_si
        
    except subprocess.TimeoutExpired:
        print("C++ program timed out")
        return None, None, None
    except json.JSONDecodeError as e:
        print(f"Failed to parse C++ output: {e}")
        print(f"Raw output: {result.stdout}")
        print(f"Stderr: {result.stderr}")
        return None, None, None
    except FileNotFoundError:
        print(f"C++ test program not found: {CPP_TEST_PROG}")
        return None, None, None


def compile_cpp():
    """Compile the C++ test program."""
    build_dir = os.path.join(TOMATOR_SRC, 'build')
    
    # Use shell=True to allow glob patterns
    compile_cmd = (
        'g++ -O2 -std=c++17 -fopenmp '
        '-I.. -I../Eigen '
        '../tests/test_single_collision.cpp '
        '../Funcs/*.cpp ../Vars/*.cpp ../Vars/*.c '
        '-o test_single_collision'
    )
    
    try:
        result = subprocess.run(
            compile_cmd,
            cwd=build_dir,
            capture_output=True,
            text=True,
            shell=True,
            timeout=120
        )
        if result.returncode != 0:
            print(f"Compilation failed:")
            print(result.stderr)
            return False
        print("C++ test program compiled successfully")
        return True
    except Exception as e:
        print(f"Compilation error: {e}")
        return False


def compare_results(py_dict, cpp_dict, name, tolerance=1e-6):
    """Compare Python and C++ results, return list of differences."""
    differences = []
    
    # Get all keys from both
    all_keys = set(py_dict.keys()) | set(cpp_dict.keys())
    
    for key in sorted(all_keys):
        py_val = py_dict.get(key, 0.0)
        cpp_val = cpp_dict.get(key, 0.0)
        
        # Skip if both are essentially zero
        if abs(py_val) < 1e-30 and abs(cpp_val) < 1e-30:
            continue
        
        # Calculate relative difference
        max_val = max(abs(py_val), abs(cpp_val))
        if max_val > 1e-30:
            rel_diff = abs(py_val - cpp_val) / max_val
        else:
            rel_diff = 0.0
        
        if rel_diff > tolerance:
            differences.append({
                'key': key,
                'py': py_val,
                'cpp': cpp_val,
                'rel_diff': rel_diff,
            })
    
    return differences


def test_single_flag(flag_name, tolerance=1e-6, verbose=True, conditions=None):
    """Test a single reaction flag, compare C++ and Python."""
    if verbose:
        print(f"\n{'='*60}")
        print(f"Testing: {flag_name}")
        print('='*60)
    
    # Run Python
    py_dn, py_dE, py_nu = run_python_collision(flag_name, conditions)
    
    # Run C++ (if available)
    cpp_dn, cpp_dE, cpp_nu = run_cpp_collision(flag_name, conditions)
    
    if cpp_dn is None:
        if verbose:
            print("  C++ not available, showing Python results only:")
            print(f"  Python dn: {format_dict(py_dn)}")
            print(f"  Python dE: {format_dict(py_dE)}")
            print(f"  Python nu: {format_dict(py_nu)}")
        return {'flag': flag_name, 'status': 'CPP_NOT_AVAILABLE', 'py_dn': py_dn, 'py_dE': py_dE, 'py_nu': py_nu}
    
    # Compare results
    dn_diffs = compare_results(py_dn, cpp_dn, 'dn', tolerance)
    dE_diffs = compare_results(py_dE, cpp_dE, 'dE', tolerance)
    nu_diffs = compare_results(py_nu, cpp_nu, 'nu', tolerance)
    
    if verbose:
        if not dn_diffs and not dE_diffs and not nu_diffs:
            print(f"  ✓ MATCH (tolerance={tolerance:.0e})")
            print(f"  Python dn: {format_dict(py_dn)}")
            print(f"  Python dE: {format_dict(py_dE)}")
            print(f"  Python nu: {format_dict(py_nu)}")
        else:
            print(f"  ✗ DIFFERENCES FOUND:")
            if dn_diffs:
                print(f"    dn differences:")
                for d in dn_diffs:
                    print(f"      {d['key']}: py={d['py']:.4e}, cpp={d['cpp']:.4e}, diff={d['rel_diff']:.2%}")
            if dE_diffs:
                print(f"    dE differences:")
                for d in dE_diffs:
                    print(f"      {d['key']}: py={d['py']:.4e}, cpp={d['cpp']:.4e}, diff={d['rel_diff']:.2%}")
            if nu_diffs:
                print(f"    nu differences:")
                for d in nu_diffs:
                    print(f"      {d['key']}: py={d['py']:.4e}, cpp={d['cpp']:.4e}, diff={d['rel_diff']:.2%}")
    
    status = 'MATCH' if not dn_diffs and not dE_diffs and not nu_diffs else 'DIFF'
    return {
        'flag': flag_name,
        'status': status,
        'dn_diffs': dn_diffs,
        'dE_diffs': dE_diffs,
        'nu_diffs': nu_diffs,
        'py_dn': py_dn,
        'py_dE': py_dE,
        'py_nu': py_nu,
        'cpp_dn': cpp_dn,
        'cpp_dE': cpp_dE,
        'cpp_nu': cpp_nu,
    }


def format_dict(d, threshold=1e-30):
    """Format a dict for printing, filtering small values."""
    filtered = {k: v for k, v in d.items() if abs(v) > threshold}
    if not filtered:
        return "{}"
    return ", ".join(f"{k}={v:.4e}" for k, v in sorted(filtered.items()))


def test_all_flags(tolerance=1e-6, conditions=None):
    """Test all reaction flags."""
    if conditions is None:
        conditions = TEST_CONDITIONS
        
    print("="*80)
    print("COMPARING C++ vs PYTHON COLLISIONS - ALL REACTION FLAGS")
    print("="*80)
    print(f"\nTest conditions:")
    for k, v in conditions.items():
        print(f"  {k} = {v}")
    print(f"\nTolerance: {tolerance:.0e}")
    
    results = []
    matches = 0
    diffs = 0
    cpp_na = 0
    
    for section, flags in ALL_FLAGS.items():
        print(f"\n{'='*60}")
        print(f"SECTION: {section}")
        print('='*60)
        
        for flag in flags:
            result = test_single_flag(flag, tolerance, verbose=True, conditions=conditions)
            results.append(result)
            
            if result['status'] == 'MATCH':
                matches += 1
            elif result['status'] == 'DIFF':
                diffs += 1
            else:
                cpp_na += 1
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Total flags tested: {len(results)}")
    print(f"  Matches: {matches}")
    print(f"  Differences: {diffs}")
    print(f"  C++ not available: {cpp_na}")
    
    if diffs > 0:
        print("\nFlags with differences:")
        for r in results:
            if r['status'] == 'DIFF':
                print(f"  - {r['flag']}")
    
    # Now run with ALL flags enabled to see cumulative effect
    print("\n" + "="*80)
    print("COMBINED TEST: ALL FLAGS ENABLED")
    print("="*80)
    all_result = test_single_flag('ALL_ENABLED', tolerance, verbose=True, conditions=conditions)
    
    return results, all_result


def load_benchmark_conditions(csv_file, target_R=91.85):
    """Load conditions from a benchmark CSV file at a specific radial position.
    
    Supports both C++ format (tmain, RadialPositions, densities in cm^-3)
    and Python format (t, R, densities in m^-3).
    """
    import pandas as pd
    
    df = pd.read_csv(csv_file)
    
    # Detect format: Python uses 't' and 'R', C++ uses 'tmain' and 'RadialPositions'
    is_python_format = 't' in df.columns and 'R' in df.columns
    
    if is_python_format:
        time_col = 't'
        r_col = 'R'
        # Python output is already in SI units (m^-3), R is in meters
        unit_factor = 1.0
        # Convert target_R from cm to m for comparison
        target_R_compare = target_R / 100.0
    else:
        time_col = 'tmain'
        r_col = 'RadialPositions'
        # C++ output is in CGS (cm^-3), R is in cm
        unit_factor = 1e6  # cm^-3 to m^-3
        target_R_compare = target_R
    
    # Get the final time step
    final_time = df[time_col].max()
    df_final = df[df[time_col] == final_time].copy()
    
    # Find R closest to target
    df_final['R_diff'] = abs(df_final[r_col] - target_R_compare)
    row = df_final.loc[df_final['R_diff'].idxmin()]
    
    print(f"Loaded benchmark conditions from {csv_file}")
    if is_python_format:
        print(f"  Format: Python (SI units)")
        print(f"  Time: {row[time_col]:.4f} s, R: {row[r_col]*100:.4f} cm ({row[r_col]:.4f} m)")
    else:
        print(f"  Format: C++ (CGS units)")
        print(f"  Time: {row[time_col]:.4f} s, R: {row[r_col]:.4f} cm")
    
    def safe_temp(T_col, default=0.5):
        """Get temperature from T column, or return default."""
        if T_col in row and row[T_col] > 0:
            return row[T_col]
        return default
    
    def safe_temp_from_E(E, n, default=0.5):
        """Calculate T = E/(n*1.5), return default if n is very small."""
        if n > 1e6:
            return E / (n * 1.5)
        return default
    
    # Extract and convert to SI (m^-3)
    # Python format already has T columns, C++ needs to compute from E/n
    if is_python_format:
        conditions = {
            'Te': safe_temp('Te', 10.0),
            'ne': row['ne'],  # Already in m^-3
            'TH': safe_temp('TH', 0.5),
            'nH': row['nH'],
            'THi': safe_temp('THi', 8.0),
            'nHi': row['nHi'],
            'TH2': safe_temp('TH2', 0.3),
            'nH2': row['nH2'],
            'TH2i': safe_temp('TH2i', 2.0),
            'nH2i': row['nH2i'],
            'TH3i': safe_temp('TH3i', 1.5),
            'nH3i': row['nH3i'],
            'THeI': safe_temp('THeI', 0.5),
            'nHeI': row['nHeI'],
            'THeII': safe_temp('THeII', 5.0),
            'nHeII': row['nHeII'],
            'THeIII': safe_temp('THeIII', 10.0),
            'nHeIII': row['nHeIII'],
            'dt': None,
            'accur': None,
        }
    else:
        conditions = {
            'Te': safe_temp_from_E(row['Ee'], row['ne'] * unit_factor, 10.0),
            'ne': row['ne'] * unit_factor,
            'TH': safe_temp_from_E(row['EH'], row['nH'] * unit_factor, 0.5),
            'nH': row['nH'] * unit_factor,
            'THi': safe_temp_from_E(row['EHi'], row['nHi'] * unit_factor, 8.0),
            'nHi': row['nHi'] * unit_factor,
            'TH2': safe_temp_from_E(row['EH2'], row['nH2'] * unit_factor, 0.3),
            'nH2': row['nH2'] * unit_factor,
            'TH2i': safe_temp_from_E(row['EH2i'], row['nH2i'] * unit_factor, 2.0),
            'nH2i': row['nH2i'] * unit_factor,
            'TH3i': safe_temp_from_E(row['EH3i'], row['nH3i'] * unit_factor, 1.5),
            'nH3i': row['nH3i'] * unit_factor,
            'THeI': safe_temp_from_E(row['EHeI'], row['nHeI'] * unit_factor, 0.5),
            'nHeI': row['nHeI'] * unit_factor,
            'THeII': safe_temp_from_E(row['EHeII'], row['nHeII'] * unit_factor, 5.0),
            'nHeII': row['nHeII'] * unit_factor,
            'THeIII': safe_temp_from_E(row['EHeIII'], row['nHeIII'] * unit_factor, 10.0),
            'nHeIII': row['nHeIII'] * unit_factor,
            'dt': None,
            'accur': None,
        }
    
    return conditions


def find_latest_benchmark_csv(directory=BENCHMARK_PYTHON_DIR):
    """Find the most recent CSV file in the benchmark directory."""
    import glob
    
    csv_files = glob.glob(os.path.join(directory, 'Res_*.csv'))
    if not csv_files:
        print(f"No CSV files found in {directory}")
        return None
    
    # Sort by modification time, most recent last
    csv_files.sort(key=os.path.getmtime)
    latest = csv_files[-1]
    print(f"Found latest benchmark CSV: {latest}")
    return latest


def main():
    parser = argparse.ArgumentParser(description='Compare C++ and Python collision implementations')
    parser.add_argument('--flag', type=str, help='Test a specific flag')
    parser.add_argument('--all', action='store_true', help='Test all flags')
    parser.add_argument('--tolerance', type=float, default=1e-6, help='Relative tolerance for comparison')
    parser.add_argument('--python-only', action='store_true', help='Run Python only (no C++ comparison)')
    parser.add_argument('--benchmark-csv', type=str, help='Load conditions from benchmark CSV file (default: latest from python benchmark)')
    parser.add_argument('--benchmark-R', type=float, default=91.85, help='Radial position to extract from benchmark (cm)')
    parser.add_argument('--no-benchmark', action='store_true', help='Use default test conditions instead of benchmark CSV')
    
    args = parser.parse_args()
    
    # Load conditions - by default use latest benchmark CSV from python folder
    conditions = None
    if args.no_benchmark:
        print("Using default test conditions (no benchmark CSV)")
    elif args.benchmark_csv:
        conditions = load_benchmark_conditions(args.benchmark_csv, args.benchmark_R)
    else:
        # Default: find and use the latest CSV from the python benchmark folder
        latest_csv = find_latest_benchmark_csv()
        if latest_csv:
            conditions = load_benchmark_conditions(latest_csv, args.benchmark_R)
    
    if args.flag:
        result = test_single_flag(args.flag, args.tolerance, conditions=conditions)
    elif args.all or True:  # Default to all
        results = test_all_flags(args.tolerance, conditions=conditions)


if __name__ == '__main__':
    main()
