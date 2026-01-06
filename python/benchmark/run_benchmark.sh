#!/bin/bash
# Benchmark script to run Python and C++ Tomator simulations
# Usage: ./run_benchmark.sh [python|cpp|both]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOMATOR_PYTHON_DIR="$SCRIPT_DIR/../tomator_dolfinx"
TOMATOR_CPP_DIR="$SCRIPT_DIR/../../src/build"

JSON_FILE="$SCRIPT_DIR/Benchmark_python_cpp.json"
PYTHON_OUTPUT_DIR="$SCRIPT_DIR/python"
CPP_OUTPUT_DIR="$SCRIPT_DIR/cpp"

# Default: run both
RUN_MODE="${1:-both}"

run_python() {
    echo "=========================================="
    echo "Running Python Tomator..."
    echo "=========================================="
    
    # Activate conda environment
    source ~/miniforge3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
    conda activate t1dl-env
    
    # Create output directory
    mkdir -p "$PYTHON_OUTPUT_DIR"
    
    # Run Python version
    cd "$TOMATOR_PYTHON_DIR"
    python examples/run_from_json.py "$JSON_FILE" -o "$PYTHON_OUTPUT_DIR"
    
    echo "Python output saved to: $PYTHON_OUTPUT_DIR"
}

run_cpp() {
    echo "=========================================="
    echo "Submitting C++ Tomator batch job..."
    echo "=========================================="
    
    # Create output directory
    mkdir -p "$CPP_OUTPUT_DIR"
    
    # Check if batch file exists
    BATCH_FILE="$CPP_OUTPUT_DIR/ITERcl_run.batch"
    if [ ! -f "$BATCH_FILE" ]; then
        echo "Error: Batch file not found at $BATCH_FILE"
        return 1
    fi
    
    # Check if C++ executable exists
    CPP_EXE="$TOMATOR_CPP_DIR/Tomator1D"
    if [ ! -f "$CPP_EXE" ]; then
        echo "Error: C++ executable not found at $CPP_EXE"
        echo "Please build the C++ code first:"
        echo "  cd $TOMATOR_CPP_DIR && cmake .. && make"
        return 1
    fi
    
    # Copy executable to cpp output directory
    cp "$CPP_EXE" "$CPP_OUTPUT_DIR/"
    
    # Submit batch job
    cd "$CPP_OUTPUT_DIR"
    sbatch ITERcl_run.batch
    
    echo "C++ batch job submitted. Output will be saved to: $CPP_OUTPUT_DIR"
}

case "$RUN_MODE" in
    python)
        run_python
        ;;
    cpp)
        run_cpp
        ;;
    both)
        run_cpp
        echo ""
        run_python
        ;;
    *)
        echo "Usage: $0 [python|cpp|both]"
        echo "  python - Run Python version only"
        echo "  cpp    - Run C++ version only"
        echo "  both   - Run both versions (default)"
        exit 1
        ;;
esac

echo ""
echo "=========================================="
echo "Benchmark complete!"
echo "=========================================="
