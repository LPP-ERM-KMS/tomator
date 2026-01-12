# Tomator-dolfinx: Python FEM Plasma Transport Solver

Python implementation of the Tomator 1D plasma transport solver using FEniCS/dolfinx.

## Overview

This package solves coupled reaction-diffusion-advection equations for plasma species (ions, neutrals, electrons) in cylindrical geometry using the finite element method.

### Governing Equation

The solver solves the cylindrical transport equation for each species:

$$\frac{\partial n}{\partial t} = \frac{1}{r} \frac{\partial}{\partial r}\left(r D \frac{\partial n}{\partial r}\right) - \frac{1}{r} \frac{\partial}{\partial r}(r V n) + S$$

where:
- $n$ = species density [m⁻³]
- $D$ = diffusion coefficient [m²/s]
- $V$ = advection velocity [m/s]
- $S$ = source/sink terms [m⁻³s⁻¹]
- $r$ = radial coordinate [m]

### Physics Steps (per timestep)

1. **BDF2 Coefficients** — Update time discretization weights
2. **Collision Rates** — Compute ionization, recombination, charge exchange rates (ADAS data)
3. **Transport Coefficients** — Calculate D (Bohm/gyro-geometric) and V for each species
4. **Parallel Losses** — Compute limiter and Bpol losses (particle sinks at boundaries)
5. **RF Power** — Add electron and ion heating from IC/EC waves
6. **Transport Solve** — Advance densities via FEM weak form
7. **Reaction Solve** — Apply atomic/molecular reactions (see solver modes below)
8. **Energy Solve** — Advance electron energy equation
9. **Post-processing** — Apply density floors, quasi-neutrality, temperature clamps

### Solver Modes

Two approaches are available for handling reaction sources:

#### 1. Explicit Mode (`operator_splitting: false`)
All sources (collisions + parallel losses) are combined and treated **explicitly** in the transport step:

$$n^{k+1} = n^k + \Delta t \left[ \nabla \cdot (D \nabla n - V n) + S^k \right]$$

- Simpler implementation
- May require smaller timesteps for stiff reactions

#### 2. Operator Splitting Mode (`operator_splitting: true`, default)
Transport and reactions are solved **separately**:

1. **Transport step**: Solve diffusion/advection for ions, neutrals, and electron energy with parallel losses only
2. **Reaction step**: Solve reaction ODE system **implicitly** at each mesh point for all densities + electron energy

$$\frac{dn}{dt} = S_{\text{reactions}}(n, T) \quad \text{(implicit Newton solve)}$$

- More robust for stiff reaction systems
- Better conservation properties
- Allows larger timesteps

### Module Overview

| Module | Description |
|--------|-------------|
| `solver.py` | BDF2 time stepper, transport equation assembly |
| `species.py` | Species definitions, PlasmaState container |
| `mesh.py` | 1D mesh generation |
| `boundary.py` | Robin/Dirichlet boundary conditions |
| `transport.py` | Diffusion models (Bohm, gyro-geometric) |
| `reactions/` | Collision rates, ADAS data, implicit reaction solver |
| `parallel.py` | Limiter and Bpol loss calculations |
| `io/` | JSON input, CSV output |
| `gui/` | Bokeh-based interactive plotter |

## Installation

### 1. conda

create and activate a conda environment:
```bash
# From the python/ folder
conda env create -f environment.yml

# Activate the environment
conda activate t1dl-env
```

### 2. pip

First install [dolfinx](https://github.com/FEniCS/dolfinx), next [create and/or
source an environment](https://docs.python.org/3/library/venv.html) and
subsequently, while in the python folder:

```bash
pip install .
```

## Verify Installation

```bash
python -c "from tomator_dolfinx import run_simulation; print('OK')"
```

## Folder Structure

```
python/
├── README.md                     # This file
├── pyproject.toml                # python project file
├── tests/                        # C++ vs Python benchmark scripts
│   └── README.md                 # Benchmark documentation
├── src/tomator_dolfinx/          # Main package
│       ├── __init__.py           # Package exports
│       ├── mesh.py               # Mesh generation (uniform, from JSON)
│       ├── species.py            # Species and PlasmaState classes
│       ├── boundary.py           # Boundary conditions
│       ├── transport.py          # Transport coefficients (Bohm, gyro-geometric)
│       ├── solver.py             # BDF2 time-stepping solver
│       ├── parallel.py           # Parallel loss calculations
│       ├── reactions/            # Collision and reaction rates
│       │   ├── collisions.py     # Collision source terms
│       │   └── rates_adas.py     # ADAS atomic data interpolation
│       ├── coupledpower/         # RF power coupling
│       ├── io/                   # Input/output handling
│       │   ├── json_input.py     # JSON input file parser
│       │   └── output.py         # CSV output writer
│       ├── gui/                  # Interactive plotting
│       └── examples/             # Example input files and scripts
│           ├── run_from_json.py  # Main entry point
│           └── *.json            # Example configurations
└── TomatorResults/               # Default output directory (git-ignored)
```

## Quick Start

### Run from JSON Input File

```bash
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
