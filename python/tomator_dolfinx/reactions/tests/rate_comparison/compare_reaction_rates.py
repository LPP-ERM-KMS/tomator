#!/usr/bin/env python3
"""
Compare C++ and Python reaction rate implementations.

This script reads C++ exported rates from CSV and compares them with
Python ReactionRates calculations, flagging discrepancies > 0.1%.

Usage: python compare_reaction_rates.py <cpp_rates.csv> [output_report.csv]
"""

import sys
import csv
import numpy as np
from pathlib import Path
from typing import Optional, Callable, Dict, Tuple
from dataclasses import dataclass

# Add the parent package to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

try:
    from tomator_dolfinx.reactions.rates import ReactionRates
    from tomator_dolfinx.reactions import rates as rates_module
except ImportError as e:
    print(f"Warning: Could not import tomator_dolfinx.reactions.rates: {e}")
    print("Some comparisons may not be available.")
    ReactionRates = None
    rates_module = None

try:
    from tomator_dolfinx.reactions import hydhel as hydhel_module
except ImportError:
    hydhel_module = None


@dataclass
class ComparisonResult:
    """Result of a single rate comparison."""
    cpp_func: str
    cpp_id: str
    param1_name: str
    param1_value: float
    param2_name: str
    param2_value: float
    cpp_rate: float
    py_rate: float
    rel_diff: float
    status: str  # 'MATCH', 'MISMATCH', 'NOT_IMPLEMENTED', 'ERROR'
    message: str


def get_python_callable(cpp_func: str, cpp_id: str) -> Optional[Callable]:
    """
    Get the Python callable for a given C++ reaction.
    
    Returns None if the reaction is not implemented in Python.
    """
    if ReactionRates is None:
        return None
    
    # Mapping from C++ IDs to Python methods
    # Format: (cpp_id) -> (method_name, needs_ne, needs_T2)
    
    # Mapping from C++ RR enum to HYDHEL reaction ID (for hydhel.compute_rate)
    # These use the generic HYDHEL polynomial fits
    # NOTE: REAC218a uses analytical formula from collisions.cpp, mapped in rr_mapping
    # NOTE: REAC218b is excluded - not used anywhere in collisions.cpp
    hydhel_id_mapping = {
        'REAC211': '2.1.1',
        'REAC212': '2.1.2',
        'REAC213': '2.1.3',   # or 2.1.5 (ionization)
        'REAC214a': '2.1.4a',
        'REAC214b': '2.1.4b',
        'REAC215': '2.1.5',
        'REAC217': '2.1.7',
        # REAC218a: Uses analytical formula, mapped in rr_mapping below
        # REAC218b: Excluded (not used in collisions.cpp)
        'REAC221a': '2.2.1a',
        'REAC221b': '2.2.1b',
        'REAC222': '2.2.2',
        'REAC223': '2.2.3',
        'REAC224': '2.2.4',
        'REAC225': '2.2.5',
        'REAC226': '2.2.6',
        'REAC227': '2.2.7',
        'REAC228': '2.2.8',
        'REAC229': '2.2.9',
        'REAC2210': '2.2.10',
        'REAC2211': '2.2.11',
        'REAC2212': '2.2.12',
        'REAC2213': '2.2.13',
        'REAC2214': '2.2.14',
        'REAC2215': '2.2.15',
        'REAC2216': '2.2.16',
    }
    
    # RR reactions (HYDHEL) - method-based where available, or use hydhel module
    rr_mapping = {
        'REAC211': ('H_excitation', False, False),  # Only Te
        'REAC212': ('H_excitation', False, False),  # Only Te (same method)
        'REAC213': ('H_ionization', False, False),
        'REAC217': ('H_radiative_recombination', False, False),
        # REAC218a: Analytical formula from collisions.cpp (same as H_radiative_recombination)
        'REAC218a': ('H_radiative_recombination', False, False),
        'REAC221a': ('H2_excitation_vib_a', False, False),
        'REAC221b': ('H2_excitation_vib_b', False, False),
        'REAC222': ('H2_excitation_elec_B', False, False),
        'REAC223': ('H2_excitation_elec_C', False, False),
        'REAC224': ('H2_excitation_elec_EF', False, False),
        'REAC2210': ('H2i_dissociative_ionization', False, False),
        'REAC2212': ('H2i_dissociation', False, False),
        'REAC2213': ('H2i_dissociation_excitation', False, False),
        'REAC2214': ('H2i_dissociative_recombination', False, False),
        'REAC2215': ('H3i_recombination_3H', False, False),
        'REAC2216': ('H3i_dissociation', False, False),
    }
    
    # For RR reactions, first check if this should use hydhel module (most RR reactions)
    if cpp_func == 'RR' and cpp_id in hydhel_id_mapping and hydhel_module is not None:
        hydhel_id = hydhel_id_mapping[cpp_id]
        # Return a lambda that calls hydhel.compute_rate with the mapped ID
        def make_hydhel_func(hid):
            def hydhel_func(Te):
                return hydhel_module.compute_rate(hid, Te)
            return hydhel_func
        return (make_hydhel_func(hydhel_id), False, False, 'hydhel')
    
    # For RR reactions NOT in hydhel_id_mapping (like REAC218a), use rr_mapping
    if cpp_func == 'RR' and cpp_id in rr_mapping:
        method_name, needs_ne, needs_T2 = rr_mapping[cpp_id]
        if hasattr(ReactionRates, method_name):
            return (getattr(ReactionRates, method_name), needs_ne, needs_T2, 'class')
    
    # RRH2 reactions (Dirk W. tables)
    rrh2_mapping = {
        'DISS': ('H2_dissociation', True, False),  # Te, ne
        'IONI': ('H2_ionization', True, False),    # Te, ne
        'RECO': ('H2i_recombination', False, False),  # Te only
        'ELAS': ('H2_elastic', False, False),       # Te only
    }
    
    # RRHe reactions (ADAS)
    rrhe_mapping = {
        'IHE1': ('HeI_ionization', True, False),   # Te, ne
        'IHE2': ('HeII_ionization', True, False),  # Te, ne
        'RHE2': ('HeII_recombination', True, False),  # Te, ne
        'RHE3': ('HeIII_recombination', True, False),  # Te, ne
    }
    
    # RRCX reactions (Charge Exchange)
    rrcx_mapping = {
        'CXHe3H': ('HeIIIH_charge_exchange', False, True),   # T1, T2
        'CXHe2H': ('HeIIH_charge_exchange', False, True),    # T1, T2
        'CXHe3He1': ('HeIIIHeI_cx_double', False, True),     # T1, T2
    }
    
    # RRion reactions (Ion-Neutral)
    rrion_mapping = {
        'REAC311': ('HiH_excitation_a', False, True),  # T1, T2
        'REAC312': ('HiH_excitation_b', False, True),
        'REAC316': ('HiH_ionization', False, True),
        'REAC321': ('HiH2_excitation_a', False, True),
        'REAC322': ('HiH2_excitation_b', False, True),
        'REAC323': ('Hi_H2_charge_exchange', False, True),
        'REAC325': ('HiH2_ionization', False, True),
        'REAC326': ('HiH2i_dissociation', False, True),
        'REAC332': ('HiHeI_ionization', False, True),
        'REAC431': ('H2iH2_charge_exchange', False, True),
        'REAC433': ('H2iH2_to_H3i', False, True),
        'REAC523': ('HeIIH2_cx_dissociation', False, True),
        'REAC531': ('HeIIHeI_charge_exchange', False, True),
        'REAC631': ('HeIIIHeI_cx_symmetric', False, True),
    }
    
    # RRel reactions (Elastic)
    # For self-collision reactions (HeIHeI, H2H2, HH), C++ only uses T1 (hard sphere formula)
    # These are marked with needs_T2=False
    rrel_mapping = {
        'HH2': ('_elastic_rate_HH2', False, True),
        'HiH': ('_elastic_rate_HiH', False, True),
        'H2iH': ('_elastic_rate_H2iH', False, True),
        'H3iH': ('_elastic_rate_H3iH', False, True),
        'HeIIH': ('_elastic_rate_HeIIH', False, True),
        'HiH2': ('_elastic_rate_HiH2', False, True),
        'H2iH2': ('_elastic_rate_H2iH2', False, True),
        'H3iH2': ('_elastic_rate_H3iH2', False, True),
        'HeIIH2': ('_elastic_rate_HeIIH2', False, True),
        'HHe': ('_elastic_rate_HHe', False, True),
        'HiHe': ('_elastic_rate_HiHe', False, True),
        'HeIIHe': ('_elastic_rate_HeIIHe', False, True),
        'HeIHeI': ('_elastic_rate_HeIHeI', False, False),  # Self-collision: only uses T1
        'H2H2': ('_elastic_rate_H2H2', False, False),       # Self-collision: only uses T1
        'HH': ('_elastic_rate_HH', False, False),           # Self-collision: only uses T1
    }
    
    # Combine all mappings
    all_mappings = {
        'RR': rr_mapping,
        'RRH2': rrh2_mapping,
        'RRHe': rrhe_mapping,
        'RRCX': rrcx_mapping,
        'RRion': rrion_mapping,
        'RRel': rrel_mapping,
    }
    
    mapping = all_mappings.get(cpp_func, {})
    if cpp_id not in mapping:
        return None
    
    method_name, needs_ne, needs_T2 = mapping[cpp_id]
    
    # Try to get the method from ReactionRates class
    if hasattr(ReactionRates, method_name):
        return (getattr(ReactionRates, method_name), needs_ne, needs_T2, 'class')
    
    # Try module-level function (for elastic rates)
    if rates_module and hasattr(rates_module, method_name):
        return (getattr(rates_module, method_name), needs_ne, needs_T2, 'module')
    
    return None


def compute_python_rate(cpp_func: str, cpp_id: str, 
                       param1_name: str, param1_value: float,
                       param2_name: str, param2_value: float) -> Tuple[float, str]:
    """
    Compute the Python rate for comparison.
    
    Returns (rate, status) where status is:
    - 'OK': Rate computed successfully
    - 'NOT_IMPLEMENTED': Method not found in Python
    - 'ERROR: <msg>': Error during computation
    """
    callable_info = get_python_callable(cpp_func, cpp_id)
    
    if callable_info is None:
        return 0.0, 'NOT_IMPLEMENTED'
    
    func, needs_ne, needs_T2, source = callable_info
    
    try:
        Te = param1_value if param1_name == 'Te' else param2_value if param2_name == 'Te' else param1_value
        T1 = param1_value if param1_name == 'T1' else param1_value
        T2_val = param2_value if param2_name == 'T2' else 1.0
        ne = param2_value if param2_name == 'ne' else 1e11
        
        # Convert scalars to arrays for methods that expect arrays
        Te_arr = np.atleast_1d(Te)
        T1_arr = np.atleast_1d(T1)
        T2_arr = np.atleast_1d(T2_val)
        ne_arr = np.atleast_1d(ne)
        
        if needs_ne:
            # 2D: Te, ne
            rate = func(Te_arr, ne_arr)
        elif needs_T2:
            # 2D: T1, T2
            rate = func(T1_arr, T2_arr)
        else:
            # 1D: Te only
            rate = func(Te_arr)
        
        # Handle array output
        if isinstance(rate, np.ndarray):
            rate = float(rate.flat[0]) if rate.size > 0 else 0.0
        
        return float(rate), 'OK'
        
    except Exception as e:
        return 0.0, f'ERROR: {str(e)}'


def compare_rates(cpp_csv_path: str, 
                  tolerance: float = 0.001,
                  output_csv: Optional[str] = None) -> Dict:
    """
    Compare C++ and Python reaction rates.
    
    Parameters
    ----------
    cpp_csv_path : str
        Path to C++ exported rates CSV file.
    tolerance : float
        Relative tolerance for matching (default 0.1% = 0.001).
    output_csv : str, optional
        Path to write detailed comparison results.
    
    Returns
    -------
    summary : dict
        Summary statistics of the comparison.
    """
    results = []
    summary = {
        'total': 0,
        'match': 0,
        'mismatch': 0,
        'not_implemented': 0,
        'error': 0,
        'by_func': {},
        'mismatches': [],
    }
    
    print(f"Reading C++ rates from: {cpp_csv_path}")
    
    with open(cpp_csv_path, 'r') as f:
        reader = csv.DictReader(f)
        
        for row in reader:
            cpp_func = row['cpp_func']
            cpp_id = row['cpp_id']
            param1_name = row['param1_name']
            param1_value = float(row['param1_value'])
            param2_name = row['param2_name']
            param2_value = float(row['param2_value']) if row['param2_value'] != '0' else 0.0
            cpp_rate = float(row['rate_cm3_s'])
            
            # Compute Python rate
            py_rate, status = compute_python_rate(
                cpp_func, cpp_id, param1_name, param1_value, param2_name, param2_value
            )
            
            # Compute relative difference
            if status == 'OK':
                if cpp_rate == 0.0 and py_rate == 0.0:
                    rel_diff = 0.0
                    result_status = 'MATCH'
                elif cpp_rate == 0.0:
                    rel_diff = float('inf') if py_rate != 0.0 else 0.0
                    result_status = 'MISMATCH' if py_rate != 0.0 else 'MATCH'
                else:
                    rel_diff = abs(py_rate - cpp_rate) / abs(cpp_rate)
                    result_status = 'MATCH' if rel_diff <= tolerance else 'MISMATCH'
            elif status == 'NOT_IMPLEMENTED':
                rel_diff = float('nan')
                result_status = 'NOT_IMPLEMENTED'
            else:
                rel_diff = float('nan')
                result_status = 'ERROR'
            
            message = status if status != 'OK' else (
                f'OK (diff={rel_diff:.2e})' if result_status == 'MATCH' 
                else f'DIFF={rel_diff:.2e}'
            )
            
            result = ComparisonResult(
                cpp_func=cpp_func,
                cpp_id=cpp_id,
                param1_name=param1_name,
                param1_value=param1_value,
                param2_name=param2_name,
                param2_value=param2_value,
                cpp_rate=cpp_rate,
                py_rate=py_rate,
                rel_diff=rel_diff,
                status=result_status,
                message=message
            )
            results.append(result)
            
            # Update summary
            summary['total'] += 1
            if result_status == 'MATCH':
                summary['match'] += 1
            elif result_status == 'MISMATCH':
                summary['mismatch'] += 1
                summary['mismatches'].append(result)
            elif result_status == 'NOT_IMPLEMENTED':
                summary['not_implemented'] += 1
            else:
                summary['error'] += 1
            
            # Per-function stats
            key = f"{cpp_func}:{cpp_id}"
            if key not in summary['by_func']:
                summary['by_func'][key] = {'match': 0, 'mismatch': 0, 'not_impl': 0, 'error': 0}
            
            if result_status == 'MATCH':
                summary['by_func'][key]['match'] += 1
            elif result_status == 'MISMATCH':
                summary['by_func'][key]['mismatch'] += 1
            elif result_status == 'NOT_IMPLEMENTED':
                summary['by_func'][key]['not_impl'] += 1
            else:
                summary['by_func'][key]['error'] += 1
    
    # Write detailed output if requested
    if output_csv:
        print(f"Writing detailed comparison to: {output_csv}")
        with open(output_csv, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([
                'cpp_func', 'cpp_id', 'param1_name', 'param1_value', 
                'param2_name', 'param2_value', 'cpp_rate', 'py_rate', 
                'rel_diff', 'status', 'message'
            ])
            for r in results:
                writer.writerow([
                    r.cpp_func, r.cpp_id, r.param1_name, r.param1_value,
                    r.param2_name, r.param2_value, r.cpp_rate, r.py_rate,
                    r.rel_diff, r.status, r.message
                ])
    
    return summary


def print_summary(summary: Dict, tolerance: float = 0.001):
    """Print a human-readable summary of the comparison."""
    
    print("\n" + "=" * 80)
    print("REACTION RATE COMPARISON SUMMARY")
    print("=" * 80)
    print(f"\nTolerance: {tolerance * 100:.2f}%")
    print(f"\nOverall Results:")
    print(f"  Total comparisons:   {summary['total']:>8d}")
    print(f"  Matches (<{tolerance*100:.1f}%):    {summary['match']:>8d}  ({100*summary['match']/max(1,summary['total']):.1f}%)")
    print(f"  Mismatches:          {summary['mismatch']:>8d}  ({100*summary['mismatch']/max(1,summary['total']):.1f}%)")
    print(f"  Not implemented:     {summary['not_implemented']:>8d}  ({100*summary['not_implemented']/max(1,summary['total']):.1f}%)")
    print(f"  Errors:              {summary['error']:>8d}  ({100*summary['error']/max(1,summary['total']):.1f}%)")
    
    print("\n" + "-" * 80)
    print("Results by Reaction:")
    print("-" * 80)
    print(f"{'Reaction':<25} {'Match':>8} {'Mismatch':>10} {'Not Impl':>10} {'Error':>8}")
    print("-" * 80)
    
    for key in sorted(summary['by_func'].keys()):
        stats = summary['by_func'][key]
        total = stats['match'] + stats['mismatch'] + stats['not_impl'] + stats['error']
        status = '✓' if stats['mismatch'] == 0 and stats['error'] == 0 else '✗' if stats['mismatch'] > 0 else '○'
        print(f"{status} {key:<23} {stats['match']:>8} {stats['mismatch']:>10} {stats['not_impl']:>10} {stats['error']:>8}")
    
    if summary['mismatches']:
        print("\n" + "-" * 80)
        print("Sample Mismatches (first 20):")
        print("-" * 80)
        for r in summary['mismatches'][:20]:
            print(f"  {r.cpp_func}:{r.cpp_id} @ {r.param1_name}={r.param1_value:.2e}")
            print(f"    C++: {r.cpp_rate:.6e}  Python: {r.py_rate:.6e}  Diff: {r.rel_diff:.2e}")
    
    print("\n" + "=" * 80)
    
    # Final verdict
    implemented = summary['match'] + summary['mismatch']
    if summary['mismatch'] == 0 and implemented > 0:
        print("✓ ALL IMPLEMENTED REACTIONS MATCH WITHIN TOLERANCE")
    elif summary['mismatch'] > 0:
        print(f"✗ {summary['mismatch']} MISMATCHES FOUND - INVESTIGATION REQUIRED")
    
    if summary['not_implemented'] > 0:
        not_impl_reactions = set()
        for key, stats in summary['by_func'].items():
            if stats['not_impl'] > 0:
                not_impl_reactions.add(key)
        print(f"\n⚠ {len(not_impl_reactions)} reactions not implemented in Python:")
        for r in sorted(not_impl_reactions):
            print(f"    - {r}")
    
    print("=" * 80)


def main():
    if len(sys.argv) < 2:
        print("Usage: python compare_reaction_rates.py <cpp_rates.csv> [output_report.csv]")
        print("\nThis script compares C++ and Python reaction rate implementations.")
        print("Run cpp_rate_export first to generate the C++ rates CSV file.")
        sys.exit(1)
    
    cpp_csv_path = sys.argv[1]
    output_csv = sys.argv[2] if len(sys.argv) > 2 else None
    
    tolerance = 0.001  # 0.1%
    
    summary = compare_rates(cpp_csv_path, tolerance=tolerance, output_csv=output_csv)
    print_summary(summary, tolerance=tolerance)
    
    # Exit with error code if mismatches found
    if summary['mismatch'] > 0:
        sys.exit(1)
    
    sys.exit(0)


if __name__ == "__main__":
    main()
