# Collision Comparison: C++ vs Python

This folder contains tools for comparing collision source terms between the C++ Tomator1D code and the Python implementation.

## Overview

The comparison validates that Python's `compute_collision_sources()` produces identical results to C++'s `collisions()` function for all reaction flags.

## Files

| File | Description |
|------|-------------|
| `compare_cpp_py_reactions.py` | Main comparison script - tests all 54 reaction flags |
| `test_collisions_cpp.cpp` | C++ test driver that calls collisions() with controlled inputs |
| `test_collisions_cpp` | Compiled C++ test binary |
| `build_test_collisions.sh` | Build script for compiling the C++ test driver |
| `CPP_VS_PYTHON_COMPARISON.md` | Detailed comparison documentation |

## Prerequisites

- Python environment with `tomator_dolfinx` package
- g++ compiler with C++17 support and OpenMP
- Access to Tomator1D source code at `../../../../../../src/`

## Building the C++ Test Driver

```bash
cd collision_comparison
./build_test_collisions.sh
```

This compiles `test_collisions_cpp.cpp` against the Tomator1D source files and produces the `test_collisions_cpp` binary.

## Running the Comparison

### Compare All Reactions (default)

```bash
python compare_cpp_py_reactions.py
```

### Compare with Custom Tolerance

```bash
# Default tolerance is 1e-6 (0.0001%)
python compare_cpp_py_reactions.py --tolerance 1e-3  # 0.1% tolerance
```

### Compare a Single Reaction Flag

```bash
python compare_cpp_py_reactions.py --flag bH_ion
python compare_cpp_py_reactions.py --flag belas_HiH
python compare_cpp_py_reactions.py --flag COULOMB_ALL
```

### Available Reaction Flags

**Atomic Hydrogen (bH)**
- `bH_exc`, `bH_ion`, `bH_3body`, `bH_rec`

**Molecular Hydrogen (bH2)**
- `bH2_elas`, `bH2_exc`, `bH2_dis`, `bH2_ion`, `bH2i_rec`, `bH2i_disrec`, `bH2_dision`, `bH3i_disrec`, `bH2i_dis`, `bH2i_disexc`, `bH3i_dis`

**Helium (bHe)**
- `bHeI_ion`, `bHeII_ion`, `bHe_cooling`, `bHeII_rec`, `bHeIII_rec`

**Charge Exchange (bcx)**
- `bcx_HiH`, `bcx_HiH2`, `bcx_H2iH2`, `bcx_HeIIH`, `bcx_HeIIHeI`, `bcx_HeIIIH`, `bcx_HeIIIHeI`

**Ion-Neutral (bion)**
- `bion_HiH_exca`, `bion_HiH_excb`, `bion_HiH_ion`, `bion_HiH2_exca`, `bion_HiH2_excb`, `bion_HiH2_325`, `bion_HiH2i_326`, `bion_H2iH2_H3i`, `bion_HiHeI_ion`, `bion_HeIIH2_cxdis`

**Elastic Collisions (belas)**
- `belas_HiH`, `belas_HH`, `belas_HH2`, `belas_HiH2`, `belas_H2iH`, `belas_H2iH2`, `belas_H3iH`, `belas_H3iH2`, `belas_H2H2`, `belas_HHeI`, `belas_HeIIH`, `belas_HiHeI`, `belas_HeIIHeI`, `belas_HeIHeI`, `belas_HeIH2`, `belas_HeIIH2`

**Coulomb Collisions (bcoulomb)**
- `COULOMB_ALL`

**Combined Test**
- `ALL_ENABLED` - Tests all reactions simultaneously

## Output Format

The comparison reports for each reaction:
- `dn`: Density source terms [m⁻³/s]
- `dE`: Energy source terms [eV/m³/s]
- `nu`: Collision frequencies [Hz]

For each species: e, H, Hi, H2, H2i, H3i, HeI, HeII, HeIII

## Example Output

```
============================================================
Testing: bH_ion
============================================================
  ✓ MATCH (tolerance=1e-06)
  Python dn: H=-2.2228e+19, Hi=2.2228e+19, e=2.2228e+19
  Python dE: H=-2.9849e+19, Hi=2.9849e+19, e=-3.0231e+20
  Python nu: H=2.8882e+03, e=7.6235e+00
```

## Test Conditions

By default, tests use conditions from the latest Python benchmark result in:
`/home/ITER/wautert/Documents/Tomator-RMA/tomator/python/benchmark/python/`

This ensures realistic plasma conditions from actual simulations.

## Troubleshooting

### Build Errors
If the build fails, ensure:
1. The Tomator source is at the expected location
2. g++ and OpenMP are available
3. Run from the `collision_comparison` directory

### Import Errors
Ensure the Python environment has access to `tomator_dolfinx`:
```bash
conda activate t1dl-env  # or your environment name
```

### Path Issues
The script uses absolute paths. If files have moved, update:
- `TOMATOR_SRC` in `compare_cpp_py_reactions.py`
- `SRC_DIR` in `build_test_collisions.sh`
