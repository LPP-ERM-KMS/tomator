# Tomator-dolfinx: Python FEM Plasma Transport Solver

Python implementation of the Tomator 1D plasma transport solver using FEniCS/dolfinx.

## Overview

This package solves coupled reaction-diffusion-advection equations for plasma species (ions, neutrals, electrons) in cylindrical geometry using the finite element method.

## Installation

### 1. Create Conda Environment

```bash
# From the python/ folder
conda env create -f environment.yml

# Activate the environment
conda activate t1dl-env
```

### 2. Verify Installation

```bash
python -c "from tomator_dolfinx import run_simulation; print('OK')"
```

## Folder Structure

```
python/
├── README.md                 # This file
├── environment.yml           # Conda environment specification
├── benchmark/                # C++ vs Python benchmark scripts
│   └── README.md             # Benchmark documentation
├── tomator_dolfinx/          # Main package
│   ├── __init__.py           # Package exports
│   ├── mesh.py               # Mesh generation (uniform, from JSON)
│   ├── species.py            # Species and PlasmaState classes
│   ├── boundary.py           # Boundary conditions
│   ├── transport.py          # Transport coefficients (Bohm, gyro-geometric)
│   ├── solver.py             # BDF2 time-stepping solver
│   ├── parallel.py           # Parallel loss calculations
│   ├── reactions/            # Collision and reaction rates
│   │   ├── collisions.py     # Collision source terms
│   │   └── rates_adas.py     # ADAS atomic data interpolation
│   ├── coupledpower/         # RF power coupling
│   ├── io/                   # Input/output handling
│   │   ├── json_input.py     # JSON input file parser
│   │   └── output.py         # CSV output writer
│   ├── gui/                  # Interactive plotting
│   └── examples/             # Example input files and scripts
│       ├── run_from_json.py  # Main entry point
│       └── *.json            # Example configurations
└── TomatorResults/           # Default output directory (git-ignored)
```

## Quick Start

### Run from JSON Input File

```bash
conda activate t1dl-env
cd tomator_dolfinx

# Run with default output directory
python examples/run_from_json.py examples/TCV5151X_fixneDV.json

# Specify custom output directory
python examples/run_from_json.py examples/TCV5151X_fixneDV.json -o ../TomatorResults/my_run

# Override simulation end time
python examples/run_from_json.py examples/TCV5151X_fixneDV.json -t 0.01

# Disable interactive plotter
python examples/run_from_json.py examples/TCV5151X_fixneDV.json --no-plot

# Verbose output
python examples/run_from_json.py examples/TCV5151X_fixneDV.json -v
```

### Use as a Library

```python
from tomator_dolfinx import (
    run_simulation,
    load_input_file,
    create_uniform_mesh,
    PlasmaState,
    BDF2Solver,
)

# Option 1: Run from input file
state = run_simulation("input.json", output_dir="results/")

# Option 2: Build simulation manually
mesh, r = create_uniform_mesh(r_min=0.6, r_max=1.0, n_elements=100)
state = PlasmaState(mesh)
state.add_helium_species()
state.initialize_from_params({"Te0": 10.0, "nHi0": 1e18})
# ... configure solver and run
```

## Input Files

The solver uses JSON input files compatible with the C++ Tomator1D format. Key sections:

| Section | Description |
|---------|-------------|
| `magnetic_field` | Bt, Bv, Bh field components |
| `toroidal_machine_geometry` | R, a, b, limiter positions |
| `neutral_pressure` | pHe, pH2 neutral pressures |
| `rf_power` | Prf, frequency, ramp settings |
| `type` | Simulation mode flags |
| `time` | Time stepping parameters |
| `mesh` | Radial mesh resolution |
| `initialvalues` | Initial density/temperature profiles |

See `examples/` folder for complete examples.

## Output

Results are saved as CSV files with columns:
- `R`: Radial position [cm]
- `t`: Time [s]
- `ne`, `nHi`, `nH2i`, `nHeI`, `nHeII`, `nHeIII`: Densities [m⁻³]
- `Te`: Electron temperature [eV]
- And many more...

## Dependencies

- **FEniCS/dolfinx** ≥0.7: Finite element framework
- **PETSc** (via petsc4py): Linear algebra backend
- **MPI** (via mpi4py): Parallel computing support
- **NumPy/SciPy**: Scientific computing
- **Bokeh** (optional): Interactive plotting
- **Pandas** (optional): Data analysis

## Comparison with C++ Version

See the [benchmark/README.md](benchmark/README.md) for instructions on running comparison benchmarks between the C++ and Python implementations.

## Troubleshooting

### Import Errors
```bash
# Ensure conda environment is activated
conda activate t1dl-env

# Verify dolfinx installation
python -c "import dolfinx; print(dolfinx.__version__)"
```

### Memory Issues
For large meshes, increase available memory or reduce mesh resolution in the input file.

### MPI Errors
The solver can run in parallel with MPI:
```bash
mpirun -n 4 python examples/run_from_json.py input.json
```
