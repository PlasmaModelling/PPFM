#!/bin/bash

# ============================================================
# Default values
# ============================================================
MAIN="main.cpp"
BUILD_DIR="build"
EXEC_DIR="executables"

# ============================================================
# Parse arguments
# ============================================================
while getopts "a:" opt; do
  case ${opt} in
    a )
      MAIN="$OPTARG"
      ;;
    \? )
      echo "Usage: $0 [-a path/to/main.cpp]"
      exit 1
      ;;
  esac
done

echo "----------------------------------------"
echo "Using main: $MAIN"
echo "----------------------------------------"

# ============================================================
# Check file exists
# ============================================================
if [ ! -f "$MAIN" ]; then
    echo "ERROR: File '$MAIN' not found"
    exit 1
fi

# ============================================================
# Clean configure (important!)
# ============================================================
echo "[1/3] Configuring..."
cmake -B $BUILD_DIR -S . -DMAIN="$MAIN" || exit 1

# ============================================================
# Build
# ============================================================
echo "[2/3] Building..."
cmake --build $BUILD_DIR -j || exit 1

# ============================================================
# Extract executable name from main
# ============================================================
FILENAME=$(basename -- "$MAIN")
NAME="${FILENAME%.*}"

EXEC_PATH="$EXEC_DIR/${NAME}.out"

# ============================================================
# Run
# ============================================================
echo "[3/3] Running: $EXEC_PATH"
echo "----------------------------------------"

if [ ! -f "$EXEC_PATH" ]; then
    echo "ERROR: Executable not found: $EXEC_PATH"
    exit 1
fi

./"$EXEC_PATH"