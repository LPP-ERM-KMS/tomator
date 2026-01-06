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
    'dt': 1e-6,
    'accur': 0.1,
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
}

def get_all_false_flags():
    """Return all individual flags set to False."""
    flags = {}
    for section_flags in ALL_FLAGS.values():
        for flag in section_flags:
            flags[flag] = False
    return flags


def run_python_collision(flag_name):
    """Run Python collision with a single flag enabled."""
    tc = TEST_CONDITIONS
    
    # All individual flags False except the one we're testing
    flags = get_all_false_flags()
    flags[flag_name] = True
    
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
        dt=tc['dt'], accur=tc['accur'],
        # Enable all section flags - individual flags control reactions
        include_H=True, include_H2=True, include_He=True,
        include_cx=True, include_ion=True, include_elastic=True,
        include_coulomb=False,  # Disable Coulomb for cleaner comparison
        **flags
    )
    
    dn, dE, nu = result
    
    # Convert to dict with scalar values
    dn_dict = {k: float(v[0]) for k, v in dn.items()}
    dE_dict = {k: float(v[0]) for k, v in dE.items()}
    nu_dict = {k: float(v[0]) for k, v in nu.items()}
    
    return dn_dict, dE_dict, nu_dict


def run_cpp_collision(flag_name):
    """
    Run C++ collision test with a single flag enabled.
    Returns dn, dE dicts parsed from C++ output.
    C++ outputs in CGS (cm^-3/s), we convert to SI (m^-3/s).
    """
    if not os.path.exists(CPP_TEST_PROG):
        print(f"C++ test program not found: {CPP_TEST_PROG}")
        print("Attempting to compile...")
        compile_cpp()
        
    if not os.path.exists(CPP_TEST_PROG):
        print("Failed to compile C++ test program")
        return None, None, None
    
    # Run C++ program with flag name as argument
    try:
        result = subprocess.run(
            [CPP_TEST_PROG, flag_name],
            capture_output=True,
            text=True,
            timeout=10
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


def test_single_flag(flag_name, tolerance=1e-6, verbose=True):
    """Test a single reaction flag, compare C++ and Python."""
    if verbose:
        print(f"\n{'='*60}")
        print(f"Testing: {flag_name}")
        print('='*60)
    
    # Run Python
    py_dn, py_dE, py_nu = run_python_collision(flag_name)
    
    # Run C++ (if available)
    cpp_dn, cpp_dE, cpp_nu = run_cpp_collision(flag_name)
    
    if cpp_dn is None:
        if verbose:
            print("  C++ not available, showing Python results only:")
            print(f"  Python dn: {format_dict(py_dn)}")
            print(f"  Python dE: {format_dict(py_dE)}")
        return {'flag': flag_name, 'status': 'CPP_NOT_AVAILABLE', 'py_dn': py_dn, 'py_dE': py_dE}
    
    # Compare results
    dn_diffs = compare_results(py_dn, cpp_dn, 'dn', tolerance)
    dE_diffs = compare_results(py_dE, cpp_dE, 'dE', tolerance)
    
    if verbose:
        if not dn_diffs and not dE_diffs:
            print(f"  ✓ MATCH (tolerance={tolerance:.0e})")
            print(f"  Python dn: {format_dict(py_dn)}")
            print(f"  Python dE: {format_dict(py_dE)}")
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
    
    status = 'MATCH' if not dn_diffs and not dE_diffs else 'DIFF'
    return {
        'flag': flag_name,
        'status': status,
        'dn_diffs': dn_diffs,
        'dE_diffs': dE_diffs,
        'py_dn': py_dn,
        'py_dE': py_dE,
        'cpp_dn': cpp_dn,
        'cpp_dE': cpp_dE,
    }


def format_dict(d, threshold=1e-30):
    """Format a dict for printing, filtering small values."""
    filtered = {k: v for k, v in d.items() if abs(v) > threshold}
    if not filtered:
        return "{}"
    return ", ".join(f"{k}={v:.4e}" for k, v in sorted(filtered.items()))


def test_all_flags(tolerance=1e-6):
    """Test all reaction flags."""
    print("="*80)
    print("COMPARING C++ vs PYTHON COLLISIONS - ALL REACTION FLAGS")
    print("="*80)
    print(f"\nTest conditions:")
    for k, v in TEST_CONDITIONS.items():
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
            result = test_single_flag(flag, tolerance, verbose=True)
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
    
    return results


def main():
    parser = argparse.ArgumentParser(description='Compare C++ and Python collision implementations')
    parser.add_argument('--flag', type=str, help='Test a specific flag')
    parser.add_argument('--all', action='store_true', help='Test all flags')
    parser.add_argument('--tolerance', type=float, default=1e-6, help='Relative tolerance for comparison')
    parser.add_argument('--python-only', action='store_true', help='Run Python only (no C++ comparison)')
    
    args = parser.parse_args()
    
    if args.flag:
        result = test_single_flag(args.flag, args.tolerance)
    elif args.all or True:  # Default to all
        results = test_all_flags(args.tolerance)


if __name__ == '__main__':
    main()
