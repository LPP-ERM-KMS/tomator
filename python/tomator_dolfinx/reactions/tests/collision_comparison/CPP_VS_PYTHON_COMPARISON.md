# C++ vs Python Collision Comparison Results

## Test Setup

Test conditions (SI units for Python, CGS for C++):
- Te = 10.0 eV, ne = 1e18 m⁻³
- TH = 0.5 eV, nH = 1e16 m⁻³
- THi = 8.0 eV, nHi = 9e17 m⁻³
- TH2 = 0.3 eV, nH2 = 5e15 m⁻³
- TH2i = 2.0 eV, nH2i = 1e15 m⁻³
- TH3i = 1.5 eV, nH3i = 1e14 m⁻³
- THeI = 0.5 eV, nHeI = 1e16 m⁻³
- THeII = 5.0 eV, nHeII = 5e15 m⁻³
- THeIII = 10.0 eV, nHeIII = 1e14 m⁻³

## Summary

| Section | Matching | Total | Match Rate |
|---------|----------|-------|------------|
| bH      | 4        | 4     | 100%       |
| bH2     | 8        | 11    | 73%        |
| bHe     | 5        | 5     | 100%*      |
| bcx     | 7        | 7     | 100%       |
| bion    | 5        | 10    | 50%        |
| belas   | 16       | 16    | 100%       |
| **Total** | **44** | **53** | **83%**   |

*bHe_cooling has 0.004% difference (numerical precision)

## Changes Made

### 2025-01 Fixes
1. **bADAS flag**: Set \`bADAS = true\` in test_single_collision.cpp to match Python default
2. **He reaction ne**: Changed from hardcoded \`ne=1.0e11\` to actual \`ne\` in C++ collisions.cpp for bHeI_ion, bHeII_ion, bHeII_rec, bHeIII_rec

## Detailed Results

### bH Section (Hydrogen atom reactions) - 100% Match ✓
- ✓ \`bH_exc\` - H excitation by electron
- ✓ \`bH_ion\` - H ionization by electron  
- ✓ \`bH_3body\` - 3-body recombination
- ✓ \`bH_rec\` - Radiative recombination

### bHe Section - 100% Match ✓
- ✓ \`bHeI_ion\` - HeI ionization (ADAS tables)
- ✓ \`bHeII_ion\` - HeII ionization (ADAS tables)
- ✓ \`bHe_cooling\` - He cooling (0.004% diff, numerical)
- ✓ \`bHeII_rec\` - HeII recombination (ADAS tables)
- ✓ \`bHeIII_rec\` - HeIII recombination (ADAS tables)

### bcx Section (Charge Exchange) - 100% Match ✓
- ✓ \`bcx_HiH\` - H+ + H charge exchange
- ✓ \`bcx_HiH2\` - H+ + H2 → H2+ + H
- ✓ \`bcx_H2iH2\` - H2+ + H2 resonant cx
- ✓ \`bcx_HeIIH\` - He+ + H → He + H+
- ✓ \`bcx_HeIIHeI\` - He+ + He resonant cx
- ✓ \`bcx_HeIIIH\` - He2+ + H → He+ + H+
- ✓ \`bcx_HeIIIHeI\` - He2+ + He → He+ + He+

### belas Section (Elastic collisions) - 100% Match ✓
- ✓ \`belas_HiH\` - H+ + H elastic
- ✓ \`belas_HH\` - H + H elastic (self-collision, zero)
- ✓ \`belas_HH2\` - H + H2 elastic
- ✓ \`belas_HiH2\` - H+ + H2 elastic
- ✓ \`belas_H2iH\` - H2+ + H elastic
- ✓ \`belas_H2iH2\` - H2+ + H2 elastic
- ✓ \`belas_H3iH\` - H3+ + H elastic
- ✓ \`belas_H3iH2\` - H3+ + H2 elastic
- ✓ \`belas_H2H2\` - H2 + H2 elastic (self-collision, zero)
- ✓ \`belas_HHeI\` - H + He elastic (same mass)
- ✓ \`belas_HeIIH\` - He+ + H elastic
- ✓ \`belas_HiHeI\` - H+ + He elastic
- ✓ \`belas_HeIIHeI\` - He+ + He elastic
- ✓ \`belas_HeIHeI\` - He + He elastic (self-collision, zero)
- ✓ \`belas_HeIH2\` - He + H2 elastic
- ✓ \`belas_HeIIH2\` - He+ + H2 elastic

### bH2 Section - 73% Match
- ✗ \`bH2_elas\` - Small difference (~0.01%) - 2D interpolation
- ✓ \`bH2_exc\` - H2 excitation
- ✓ \`bH2_dis\` - H2 dissociation
- ✗ \`bH2_ion\` - Small difference (~0.19%) - 2D interpolation
- ✓ \`bH2i_rec\` - H2+ recombination
- ✓ \`bH2i_disrec\` - H2+ dissociative recombination
- ✓ \`bH2_dision\` - H2 dissociative ionization
- ✓ \`bH3i_disrec\` - H3+ dissociative recombination
- ✗ \`bH2i_dis\` - 10% difference (needs investigation)
- ✓ \`bH2i_disexc\` - H2+ dissociative excitation
- ✓ \`bH3i_dis\` - H3+ dissociation

### bion Section - 50% Match
- ✓ \`bion_HiH_exca\` - H+ + H excitation
- ✓ \`bion_HiH_excb\` - H+ + H excitation (zero)
- ✓ \`bion_HiH_ion\` - H+ + H ionization (zero)
- ✓ \`bion_HiH2_exca\` - H+ + H2 excitation
- ✓ \`bion_HiH2_excb\` - H+ + H2 excitation
- ✗ \`bion_HiH2_325\` - Different dE attribution (electron vs Hi)
- ✗ \`bion_HiH2i_326\` - Different dn/dE
- ✗ \`bion_H2iH2_H3i\` - 75% difference
- ✗ \`bion_HiHeI_ion\` - Different dE attribution
- ✗ \`bion_HeIIH2_cxdis\` - Various differences

## Running the Comparison

\`\`\`bash
cd /path/to/tomator/python/tomator_dolfinx/reactions/tests
conda activate t1dl-env
python compare_cpp_py_reactions.py --all
\`\`\`

For a single flag:
\`\`\`bash
python compare_cpp_py_reactions.py --flag bH_ion
\`\`\`
