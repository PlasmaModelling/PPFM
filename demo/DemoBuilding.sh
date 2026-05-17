#!/usr/bin/env bash

# bash script for generating all demo executables.

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$( dirname "$SCRIPT_DIR" )"

cd "$PROJECT_ROOT"

for demo_file in demo/*.cpp; do
    demo_name=$(basename "$demo_file" .cpp)
    rm -rf build
    cmake -B build -S . -DMAIN="demo/${demo_name}.cpp"
    cmake --build build
done
