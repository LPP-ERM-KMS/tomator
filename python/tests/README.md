# Tomator Benchmark: C++ vs Python

This folder contains scripts and configurations to benchmark the C++ and Python implementations of Tomator1D.

## Folder Structure

```
benchmark/
├── README.md                    # This file
├── Benchmark_python_cpp.json    # Input configuration for benchmarking
├── Benchmark_cpp_python.json    # Alternative input configuration
├── run_benchmark.sh             # Main script to run benchmarks
├── plot_results.sh              # Script to visualize results
├── cpp/                         # C++ output directory
│   ├── ITERcl_run.batch         # SLURM batch file for C++ runs
│   ├── Benchmark_cpp_python/    # C++ simulation results
│   └── Log_benchmark.log        # C++ run log
└── python/                      # Python output directory
    └── Res_*.csv                # Python simulation results
```

## Prerequisites

### C++ Version
1. Build the C++ Tomator1D executable:
   ```bash
   cd ../../src/build
   cmake ..
   make
   ```

### Python Version
1. Activate the conda environment:
   ```bash
   conda activate t1dl-env
   ```

2. Ensure all dependencies are installed (see `python/environment.yml`).

## Running the Benchmark

### Quick Start
```bash
# Run both C++ and Python benchmarks
./run_benchmark.sh both

# Run Python only
./run_benchmark.sh python

# Run C++ only (submits SLURM job)
./run_benchmark.sh cpp
```

### Manual Execution

#### Python
```bash
conda activate t1dl-env
cd ../tomator_dolfinx
python -u examples/run_from_json.py ../benchmark/Benchmark_python_cpp.json -o ../benchmark/python
```

#### C++ (Local)
```bash
cd cpp
cp ../../src/build/Tomator1D .
./Tomator1D ../Benchmark_cpp_python.json
```

#### C++ (SLURM Cluster)
```bash
cd cpp
cp ../../src/build/Tomator1D .
sbatch ITERcl_run.batch
```

## Plotting Results

Use the `plot_results.sh` script to visualize simulation results:

```bash
# Plot both C++ and Python results
./plot_results.sh both

# Plot latest C++ result only
./plot_results.sh cpp

# Plot latest Python result only
./plot_results.sh python
```

This opens an interactive web-based plotter in your browser.

## Configuration Files

- **Benchmark_python_cpp.json**: Main benchmark configuration (TCV-like parameters)
  - Fixed power fraction mode (`bfixpowerfrac: true`)
  - He and H2 neutral pressures
  - IC heating at 82.7 MHz

- **Benchmark_cpp_python.json**: Alternative configuration for C++ runs

## Output Files

- **Python**: Results are saved as `Res_YYYYMMDD_HHMMSS.csv` in the `python/` folder
- **C++**: Results are saved in `cpp/Benchmark_cpp_python/` folder

## Notes

- The C++ version uses SLURM for job submission on clusters
- Results CSV files are ignored by git (see `.gitignore`)
- The `Tomator1D` binary is also ignored by git
