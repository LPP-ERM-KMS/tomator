# Reaction Rate Comparison Tests

This folder contains scripts to compare the Python reaction rate implementations
against the reference C++ implementations in `tomator/src/Vars/reactionrates.cpp`.

## Overview

The comparison verifies that all Python reaction rate functions produce results
within 0.1% of the C++ implementations across a comprehensive grid of temperatures
and densities.

## Quick Start

### Option 1: Run the Shell Script (Recommended)

```bash
cd tomator/python/tomator_dolfinx/reactions/tests
./run_rate_comparison.sh
```

This will automatically:
1. Build the C++ rate export driver
2. Export C++ rates to CSV
3. Run Python comparison
4. Generate a detailed report

### Option 2: Run Steps Manually

**Step 1: Build the C++ rate export driver**
```bash
./build_cpp_rate_export.sh
```

**Step 2: Export C++ rates to CSV**
```bash
./cpp_rate_export ../hydhel.tex cpp_rates.csv
```

**Step 3: Run Python comparison**
```bash
# Make sure you're in the t1dl-env conda environment
conda activate t1dl-env
python compare_reaction_rates.py cpp_rates.csv comparison_report.csv
```

## Files

| File | Description |
|------|-------------|
| `compare_reaction_rates.py` | Main Python comparison script |
| `cpp_rate_export.cpp` | C++ driver to export rates to CSV |
| `build_cpp_rate_export.sh` | Script to build the C++ driver |
| `run_rate_comparison.sh` | All-in-one script to run the full comparison |
| `../hydhel.tex` | HYDHEL polynomial fit coefficients (in parent folder) |

## Output Files

After running the comparison:

| File | Description |
|------|-------------|
| `cpp_rates.csv` | Exported C++ reaction rates (13,950 rows) |
| `comparison_report.csv` | Detailed comparison with C++/Python values and relative differences |

The comparison report contains columns:
- `cpp_func`, `cpp_id`: Reaction function and identifier
- `param1_name/value`, `param2_name/value`: Input parameters (Te, ne, T1, T2)
- `cpp_rate`, `py_rate`: Computed rates from each implementation
- `rel_diff`: Relative difference between implementations
- `status`: MATCH, MISMATCH, NOT_IMPLEMENTED, or ERROR
- `message`: Additional details

## Reactions Tested

The comparison covers the following reaction categories:

### RR - HYDHEL Polynomial Fits (electron-impact)
- H excitation, ionization, recombination
- H₂ vibrational/electronic excitation
- H₂⁺ dissociation and recombination
- H₃⁺ reactions

### RRH2 - Dirk Wünderlich Tables
- H₂ dissociation, ionization, recombination
- H₂ elastic collisions

### RRHe - ADAS Tables
- He⁺/He²⁺ ionization
- He²⁺/He³⁺ recombination

### RRCX - Charge Exchange
- He²⁺ + H → He⁺ + H⁺
- He³⁺ + H → He²⁺ + H⁺
- He³⁺ + He → He²⁺ + He⁺

### RRion - Ion-Neutral Reactions
- H⁺ + H/H₂ excitation, ionization
- H₂⁺ + H₂ charge exchange
- He²⁺/He³⁺ + H₂/He reactions

### RRel - Elastic Collisions
- Langevin and hard-sphere collision rates
- Ion-neutral momentum transfer

## Expected Output

When all rates match within tolerance:

```
================================================================================
✓ ALL IMPLEMENTED REACTIONS MATCH WITHIN TOLERANCE
================================================================================
```

## Requirements

- C++ compiler (g++ with C++17 support)
- Python 3.8+ with NumPy
- conda environment `t1dl-env` with tomator_dolfinx installed

## Troubleshooting

### HYDHEL file not found
The `hydhel.tex` file should be in the parent `reactions/` folder. If missing, copy it from:
```bash
cp /path/to/tomator/src/SimParams/Public/hydhel.tex ../
```

Or pass a custom path as argument:
```bash
./run_rate_comparison.sh /path/to/hydhel.tex
```

### Python import errors
Make sure the conda environment is activated:
```bash
conda activate t1dl-env
```

And that PYTHONPATH includes the tomator_dolfinx package:
```bash
export PYTHONPATH=/path/to/tomator/python:$PYTHONPATH
```
