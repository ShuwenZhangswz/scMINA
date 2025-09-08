#!/bin/bash
# Simple setup script for scMINA conda environment

set -e

echo "scMINA Setup"
echo "==============="
echo "Setting up scMINA with Python and R support"
echo ""

# Check if we're in the right directory
if [ ! -f "environment.yml" ]; then
    echo "Please run this script from the scMINA root directory"
    exit 1
fi

echo "📁 Current directory: $(pwd)"
echo ""

echo "🔧 Setting up conda environment..."
./setup_conda_env.sh

echo ""
echo "Setup complete!"
echo ""
echo "Next steps:"
echo "Activate environment: conda activate scmina"
