#!/usr/bin/env bash

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Running all .out files in: $SCRIPT_DIR"
echo

for exe in "$SCRIPT_DIR"/*.out; do
    # Handle case with no .out files
    if [[ ! -e "$exe" ]]; then
        echo "No .out files found."
        exit 0
    fi

    if [[ -x "$exe" ]]; then
        echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
        echo "Running: $exe"
        echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
        "$exe"
        echo
    else
        echo "Skipping (not executable): $exe"
    fi
done

echo "Done."
