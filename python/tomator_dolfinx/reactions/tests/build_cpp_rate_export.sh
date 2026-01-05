#!/bin/bash
#
# Build script for the C++ rate export driver
#
# This script compiles the cpp_rate_export.cpp program which links against
# the Tomator C++ source to export reaction rates for comparison with Python.
#
# Prerequisites:
# - C++ compiler (g++ or clang++)
# - Tomator source code in ../../../src/
#

set -e  # Exit on error

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
TOMATOR_SRC="$SCRIPT_DIR/../../../../src"

# Check if Tomator source exists
if [ ! -d "$TOMATOR_SRC" ]; then
    echo "Error: Tomator source directory not found at: $TOMATOR_SRC"
    exit 1
fi

echo "Building cpp_rate_export..."
echo "  Source dir: $TOMATOR_SRC"
echo "  Output:     $SCRIPT_DIR/cpp_rate_export"

# Compile the export driver
# We need to link against the relevant Tomator source files
g++ -O2 -std=c++17 \
    -I"$TOMATOR_SRC" \
    -I"$TOMATOR_SRC/SimParams/Public" \
    -I"$TOMATOR_SRC/Funcs" \
    -I"$TOMATOR_SRC/Vars" \
    "$SCRIPT_DIR/cpp_rate_export.cpp" \
    "$TOMATOR_SRC/Vars/reactionrates.cpp" \
    "$TOMATOR_SRC/Vars/simparam.cpp" \
    "$TOMATOR_SRC/Vars/constants.c" \
    -o "$SCRIPT_DIR/cpp_rate_export" \
    -lm

echo "Build successful: $SCRIPT_DIR/cpp_rate_export"
