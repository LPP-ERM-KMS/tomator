#!/bin/bash
#
# Run the complete C++ vs Python reaction rate comparison
#
# This script:
# 1. Builds the C++ rate export driver
# 2. Runs the C++ export to generate rates CSV
# 3. Runs the Python comparison script
#
# Usage: ./run_rate_comparison.sh [hydhel_path]
#
# If hydhel_path is not specified, uses the default location.
#

set -e  # Exit on error

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Default path to hydhel.tex (in parent reactions folder)
DEFAULT_HYDHEL="$SCRIPT_DIR/../hydhel.tex"
HYDHEL_PATH="${1:-$DEFAULT_HYDHEL}"

# Output paths
CPP_EXPORT="$SCRIPT_DIR/cpp_rate_export"
CSV_OUTPUT="$SCRIPT_DIR/cpp_rates.csv"
REPORT_OUTPUT="$SCRIPT_DIR/comparison_report.csv"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "=========================================="
echo "  C++ vs Python Reaction Rate Comparison"
echo "=========================================="
echo ""

# Check if hydhel.tex exists
if [ ! -f "$HYDHEL_PATH" ]; then
    echo -e "${RED}Error: hydhel.tex not found at: $HYDHEL_PATH${NC}"
    echo "Please provide the path as an argument: ./run_rate_comparison.sh /path/to/hydhel.tex"
    exit 1
fi

echo -e "${YELLOW}Step 1: Building C++ rate export driver...${NC}"
if [ ! -x "$CPP_EXPORT" ] || [ "$SCRIPT_DIR/cpp_rate_export.cpp" -nt "$CPP_EXPORT" ]; then
    bash "$SCRIPT_DIR/build_cpp_rate_export.sh"
else
    echo "  (Using existing binary - up to date)"
fi
echo ""

echo -e "${YELLOW}Step 2: Exporting C++ reaction rates...${NC}"
echo "  HYDHEL path: $HYDHEL_PATH"
echo "  Output CSV:  $CSV_OUTPUT"
"$CPP_EXPORT" "$HYDHEL_PATH" "$CSV_OUTPUT"
echo ""

echo -e "${YELLOW}Step 3: Running Python comparison...${NC}"
python3 "$SCRIPT_DIR/compare_reaction_rates.py" "$CSV_OUTPUT" "$REPORT_OUTPUT"

echo ""
echo -e "${GREEN}Comparison complete!${NC}"
echo "  Detailed report: $REPORT_OUTPUT"
