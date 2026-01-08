#!/bin/bash
# =============================================================================
# Plot Results Script for Tomator Benchmark
# =============================================================================
# Usage:
#   ./plot_results.sh cpp      - Plot latest C++ result
#   ./plot_results.sh python   - Plot latest Python result
#   ./plot_results.sh both     - Plot both (opens two browser tabs)
#   ./plot_results.sh          - Same as 'both'
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CPP_OUTPUT_DIR="$SCRIPT_DIR/cpp"
PYTHON_OUTPUT_DIR="$SCRIPT_DIR/python"

# Conda environment
CONDA_ENV="t1dl-env"

# Find an available port (random in range 5100-5999)
find_available_port() {
    local port
    for _ in {1..50}; do
        port=$((5100 + RANDOM % 900))
        # Check if port is available using /dev/tcp
        if ! (echo >/dev/tcp/localhost/$port) 2>/dev/null; then
            echo "$port"
            return 0
        fi
    done
    # Fallback: return a random port and hope for the best
    echo "$((5100 + RANDOM % 900))"
}

# Find the latest Res_*.csv file in a directory (searches subdirectories too)
find_latest_result() {
    local dir="$1"
    if [ ! -d "$dir" ]; then
        echo ""
        return
    fi
    # Find all Res_*.csv files, sort by modification time, get the newest
    find "$dir" -name "Res_*.csv" -type f -printf '%T@ %p\n' 2>/dev/null | \
        sort -n | tail -1 | cut -d' ' -f2-
}

# Launch plotter for a result file
launch_plotter() {
    local csv_file="$1"
    local label="$2"
    local port="$3"
    
    if [ -z "$csv_file" ] || [ ! -f "$csv_file" ]; then
        echo "Error: No result file found for $label"
        return 1
    fi
    
    echo "Launching plotter for $label: $csv_file"
    echo "  URL: http://localhost:$port"
    
    # Activate conda and run plotter
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate "$CONDA_ENV"
    
    cd "$SCRIPT_DIR/.."
    TOMATOR_CSV_FILE="$csv_file" python -m bokeh serve --show \
        --port "$port" \
        tomator_dolfinx/gui/plot_app.py &
}

# Main
MODE="${1:-both}"

case "$MODE" in
    cpp)
        CPP_FILE=$(find_latest_result "$CPP_OUTPUT_DIR")
        if [ -z "$CPP_FILE" ]; then
            echo "No C++ result files found in $CPP_OUTPUT_DIR"
            exit 1
        fi
        PORT=$(find_available_port)
        launch_plotter "$CPP_FILE" "C++" "$PORT"
        ;;
    
    python)
        PYTHON_FILE=$(find_latest_result "$PYTHON_OUTPUT_DIR")
        if [ -z "$PYTHON_FILE" ]; then
            echo "No Python result files found in $PYTHON_OUTPUT_DIR"
            exit 1
        fi
        PORT=$(find_available_port)
        launch_plotter "$PYTHON_FILE" "Python" "$PORT"
        ;;
    
    both)
        CPP_FILE=$(find_latest_result "$CPP_OUTPUT_DIR")
        PYTHON_FILE=$(find_latest_result "$PYTHON_OUTPUT_DIR")
        
        if [ -z "$CPP_FILE" ] && [ -z "$PYTHON_FILE" ]; then
            echo "No result files found in either directory"
            exit 1
        fi
        
        if [ -n "$CPP_FILE" ]; then
            PORT=$(find_available_port)
            launch_plotter "$CPP_FILE" "C++" "$PORT"
        else
            echo "No C++ result files found"
        fi
        
        if [ -n "$PYTHON_FILE" ]; then
            PORT=$(find_available_port)
            launch_plotter "$PYTHON_FILE" "Python" "$PORT"
        else
            echo "No Python result files found"
        fi
        ;;
    
    *)
        echo "Usage: $0 [cpp|python|both]"
        echo "  cpp    - Plot latest C++ result"
        echo "  python - Plot latest Python result"
        echo "  both   - Plot both results"
        exit 1
        ;;
esac

echo ""
echo "Plotter(s) running in background. Press Ctrl+C to stop all."
wait
