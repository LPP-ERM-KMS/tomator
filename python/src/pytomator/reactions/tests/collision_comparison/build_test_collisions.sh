#!/bin/bash
# Build script for test_collisions_cpp
# This compiles the test driver that calls the real collisions() function

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="${SCRIPT_DIR}/../../../../../src"
BUILD_DIR="${SCRIPT_DIR}"

echo "Building test_collisions_cpp..."
echo "  Source dir: $SRC_DIR"
echo "  Build dir: $BUILD_DIR"

# Check if source dir exists
if [ ! -d "$SRC_DIR" ]; then
    echo "ERROR: Source directory not found: $SRC_DIR"
    exit 1
fi

# Compile - link against all necessary Tomator source files
g++ -O2 -std=c++17 -fopenmp \
    -I"${SRC_DIR}" \
    -I"${SRC_DIR}/Eigen" \
    -I"${SRC_DIR}/Funcs" \
    -I"${SRC_DIR}/Vars" \
    "${SCRIPT_DIR}/test_collisions_cpp.cpp" \
    "${SRC_DIR}/Funcs/collisions.cpp" \
    "${SRC_DIR}/Funcs/functions.cpp" \
    "${SRC_DIR}/Vars/reactionrates.cpp" \
    "${SRC_DIR}/Vars/simparam.cpp" \
    "${SRC_DIR}/Vars/globalVariables.cpp" \
    "${SRC_DIR}/Vars/positions.cpp" \
    "${SRC_DIR}/Vars/constants.c" \
    -o "${BUILD_DIR}/test_collisions_cpp" \
    -lm

echo "Build successful!"
echo "Executable: ${BUILD_DIR}/test_collisions_cpp"
echo ""
echo "Usage:"
echo "  ${BUILD_DIR}/test_collisions_cpp <hydhel.tex path> [Te] [ne_cgs]"
echo ""
echo "Example:"
echo "  ${BUILD_DIR}/test_collisions_cpp /path/to/HYDHEL.TEX 10 1e12"
