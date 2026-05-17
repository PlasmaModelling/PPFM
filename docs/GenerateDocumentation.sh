#!/usr/bin/env bash

set -e

echo "Generating PPFM Doxygen documentation..."

doxygen Doxyfile.in

echo ""
echo "Documentation generated at:"
echo "docs/doxygen/html/index.html"
