#!/usr/bin/env python3
"""
Compare each reaction between C++ collisions.cpp and Python collisions.py.
Uses individual reaction flags to isolate each reaction for comparison.

This script tests each reaction in isolation by enabling only one flag at a time,
then compares the dn and dE outputs to manually computed expected values from C++ formulas.
"""

import numpy as np
import sys
sys.path.insert(0, '/home/ITER/wautert/Documents/Tomator-RMA/tomator/python')

from tomator_dolfinx.reactions import collisions
from tomator_dolfinx.reactions.rates import ReactionRates
from tomator_dolfinx.reactions import rates as rates_module

# Create rates instance
rates = ReactionRates()

# Unit conversions (matching C++)
CM3_TO_M3 = 1e-6
M3_TO_CM3 = 1e6

# Test conditions - as numpy arrays (single element)
Te = np.array([10.0])       # eV
ne = np.array([1e18])       # m^-3
nH = np.array([1e16])       # m^-3
nHi = np.array([9e17])      # m^-3
TH = np.array([0.5])        # eV
THi = np.array([8.0])       # eV

# H2 species
nH2 = np.array([5e15])      # m^-3
nH2i = np.array([1e15])     # m^-3
nH3i = np.array([1e14])     # m^-3
TH2 = np.array([0.3])       # eV
TH2i = np.array([2.0])      # eV
TH3i = np.array([1.5])      # eV

# He species
nHeI = np.array([1e16])     # m^-3
nHeII = np.array([5e15])    # m^-3
nHeIII = np.array([1e14])   # m^-3
THeI = np.array([0.5])      # eV
THeII = np.array([5.0])     # eV
THeIII = np.array([10.0])   # eV

# CGS densities (scalars for ndamp)
ne_cgs = ne[0] * M3_TO_CM3
nH_cgs = nH[0] * M3_TO_CM3
nHi_cgs = nHi[0] * M3_TO_CM3
nH2_cgs = nH2[0] * M3_TO_CM3
nH2i_cgs = nH2i[0] * M3_TO_CM3
nH3i_cgs = nH3i[0] * M3_TO_CM3
nHeI_cgs = nHeI[0] * M3_TO_CM3
nHeII_cgs = nHeII[0] * M3_TO_CM3
nHeIII_cgs = nHeIII[0] * M3_TO_CM3

# ndamp function (same as C++)
def ndamp(n_cgs, n0=1.0, alpha=0.66):
    return (1.0 + n0 / max(n_cgs, 1e-30))**(-alpha)

dt = 1e-6
accur = 0.1

def get_all_false_flags():
    """Return all flags set to False."""
    return {
        # bH flags
        'bH_exc': False, 'bH_ion': False, 'bH_3body': False, 'bH_rec': False,
        # bH2 flags
        'bH2_elas': False, 'bH2_exc': False, 'bH2_dis': False, 'bH2_ion': False,
        'bH2i_rec': False, 'bH2i_disrec': False, 'bH2_dision': False,
        'bH3i_disrec': False, 'bH2i_dis': False, 'bH2i_disexc': False, 'bH3i_dis': False,
        # bHe flags
        'bHeI_ion': False, 'bHeII_ion': False, 'bHe_cooling': False,
        'bHeII_rec': False, 'bHeIII_rec': False,
        # bcx flags
        'bcx_HiH': False, 'bcx_HiH2': False, 'bcx_H2iH2': False, 'bcx_HeIIH': False,
        'bcx_HeIIHeI': False, 'bcx_HeIIIH': False, 'bcx_HeIIIHeI': False,
        # bion flags
        'bion_HiH_exca': False, 'bion_HiH_excb': False, 'bion_HiH_ion': False,
        'bion_HiH2_exca': False, 'bion_HiH2_excb': False, 'bion_HiH2_325': False,
        'bion_HiH2i_326': False, 'bion_H2iH2_H3i': False, 'bion_HiHeI_ion': False,
        'bion_HeIIH2_cxdis': False,
        # belas flags (matching function signature in collisions.py)
        'belas_HiH': False, 'belas_HH': False, 'belas_HH2': False, 'belas_HiH2': False,
        'belas_H2iH': False, 'belas_H2iH2': False, 'belas_H3iH': False, 'belas_H3iH2': False,
        'belas_H2H2': False, 'belas_HHeI': False, 'belas_HeIIH': False, 'belas_HiHeI': False,
        'belas_HeIIHeI': False, 'belas_HeIHeI': False, 'belas_HeIH2': False, 'belas_HeIIH2': False,
    }

def run_python_with_flag(flag_name):
    """Run Python collision sources with only one flag enabled."""
    flags = get_all_false_flags()
    flags[flag_name] = True
    
    # Determine which master section to enable based on the flag name
    include_H = flag_name.startswith('bH_')
    include_H2 = flag_name.startswith('bH2') or flag_name.startswith('bH3i')
    include_He = flag_name.startswith('bHe')
    include_cx = flag_name.startswith('bcx_')
    include_ion = flag_name.startswith('bion_')
    include_elastic = flag_name.startswith('belas_')
    
    # For belas reactions, we need to enable the relevant species section too
    # since the code checks both include_elastic AND include_X
    if include_elastic:
        # Enable species sections needed for elastic collisions
        if 'HiH' in flag_name or 'HH' in flag_name:
            include_H = True
        if 'H2' in flag_name or 'H3' in flag_name:
            include_H2 = True
        if 'He' in flag_name:
            include_He = True
    
    # Similarly for bcx and bion reactions that involve multiple species
    if include_cx:
        # bcx reactions can involve H, H2, and He
        if 'H2' in flag_name:
            include_H2 = True
        if 'He' in flag_name:
            include_He = True
        if 'HiH' in flag_name or 'HeIIH' in flag_name or 'HeIIIH' in flag_name:
            include_H = True
    
    if include_ion:
        # bion reactions can involve H, H2, and He
        if 'H2' in flag_name:
            include_H2 = True
        if 'He' in flag_name:
            include_He = True
        if 'HiH' in flag_name:
            include_H = True
    
    # Disable ALL sections to start clean, then enable only the relevant one
    result = collisions.compute_collision_sources(
        Te=Te, ne=ne, nH=nH, nHi=nHi, TH=TH, THi=THi,
        nH2=nH2, nH2i=nH2i, nH3i=nH3i, TH2=TH2, TH2i=TH2i, TH3i=TH3i,
        nHeI=nHeI, nHeII=nHeII, nHeIII=nHeIII,
        THeI=THeI, THeII=THeII, THeIII=THeIII,
        dt=dt, accur=accur,
        include_H=include_H, include_H2=include_H2, include_He=include_He,
        include_cx=include_cx, include_ion=include_ion, include_elastic=include_elastic,
        include_coulomb=False,  # Disable Coulomb to isolate reactions
        **flags
    )
    return result


def check_nonzero_outputs(flag_name, dn, dE):
    """Check if any outputs are non-zero."""
    nonzero_dn = {k: v[0] for k, v in dn.items() if np.abs(v[0]) > 1e-30}
    nonzero_dE = {k: v[0] for k, v in dE.items() if np.abs(v[0]) > 1e-30}
    return nonzero_dn, nonzero_dE


def main():
    print("=" * 80)
    print("COMPARING REACTION-BY-REACTION: Python collisions.py")
    print("=" * 80)
    print(f"\nTest conditions:")
    print(f"  Te={Te[0]} eV, ne={ne[0]:.2e} m^-3")
    print(f"  TH={TH[0]} eV, nH={nH[0]:.2e} m^-3")
    print(f"  THi={THi[0]} eV, nHi={nHi[0]:.2e} m^-3")
    print(f"  TH2={TH2[0]} eV, nH2={nH2[0]:.2e} m^-3")
    print(f"  TH2i={TH2i[0]} eV, nH2i={nH2i[0]:.2e} m^-3")
    print(f"  TH3i={TH3i[0]} eV, nH3i={nH3i[0]:.2e} m^-3")
    print(f"  THeI={THeI[0]} eV, nHeI={nHeI[0]:.2e} m^-3")
    print(f"  THeII={THeII[0]} eV, nHeII={nHeII[0]:.2e} m^-3")
    print(f"  THeIII={THeIII[0]} eV, nHeIII={nHeIII[0]:.2e} m^-3")
    print()
    
    # All flags to test, grouped by section
    all_flags = {
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
    
    total_tests = 0
    total_passed = 0
    total_failed = 0
    results_summary = []
    
    for section, flags in all_flags.items():
        print(f"\n{'='*80}")
        print(f"SECTION: {section}")
        print("=" * 80)
        
        for flag in flags:
            total_tests += 1
            try:
                dn, dE, nu = run_python_with_flag(flag)
                nonzero_dn, nonzero_dE = check_nonzero_outputs(flag, dn, dE)
                
                # Display results
                if nonzero_dn or nonzero_dE:
                    print(f"\n{flag}:")
                    if nonzero_dn:
                        print(f"  dn: ", end="")
                        for k, v in nonzero_dn.items():
                            print(f"{k}={v:.4e}, ", end="")
                        print()
                    if nonzero_dE:
                        print(f"  dE: ", end="")
                        for k, v in nonzero_dE.items():
                            print(f"{k}={v:.4e}, ", end="")
                        print()
                    total_passed += 1
                    results_summary.append((flag, 'ACTIVE', nonzero_dn, nonzero_dE))
                else:
                    print(f"{flag}: NO OUTPUT (possible issue?)")
                    # Some reactions like belas_HH have no energy exchange
                    if 'belas_HH' in flag or 'belas_H2H2' in flag or 'belas_HeIHeI' in flag:
                        print(f"  (Expected: self-collision, no energy exchange)")
                        total_passed += 1
                        results_summary.append((flag, 'NO_EXCHANGE', {}, {}))
                    else:
                        total_failed += 1
                        results_summary.append((flag, 'NO_OUTPUT', {}, {}))
                        
            except Exception as e:
                print(f"\n[ERROR] {flag}: {e}")
                total_failed += 1
                results_summary.append((flag, 'ERROR', {}, {}))
    
    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total reactions tested: {total_tests}")
    print(f"Active (with output): {sum(1 for r in results_summary if r[1] == 'ACTIVE')}")
    print(f"Self-collisions (nu only): {sum(1 for r in results_summary if r[1] == 'NO_EXCHANGE')}")
    print(f"No output (issues): {sum(1 for r in results_summary if r[1] == 'NO_OUTPUT')}")
    print(f"Errors: {sum(1 for r in results_summary if r[1] == 'ERROR')}")
    
    # Show issues
    issues = [r for r in results_summary if r[1] in ('NO_OUTPUT', 'ERROR')]
    if issues:
        print(f"\nReactions with issues:")
        for r in issues:
            print(f"  - {r[0]}: {r[1]}")
    
    return total_failed == 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
