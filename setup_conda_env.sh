#!/bin/bash
# Setup script for scMINA conda environment with both Python and R

set -e

echo "🔧 Setting up scMINA conda environment (Python + R)"
echo "=================================================="

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "Installing miniforge..."
    cd $HOME
    curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
    bash Miniforge3-$(uname)-$(uname -m).sh -b -p $HOME/miniforge3
    source $HOME/miniforge3/bin/activate
    echo "Miniforge installed and activated"
fi

# Always use conda to avoid conflicting mamba from other toolchains
CONDA_CMD="conda"
echo "Using conda"

# Create the environment
echo "Creating scmina conda environment..."
$CONDA_CMD env create -f environment.yml

echo ""
echo "✅ scmina environment created successfully!"
echo ""
echo "To activate the environment:"
echo "  conda activate scmina"
echo ""
echo "To deactivate:"
echo "  conda deactivate"
echo ""
echo "To remove the environment:"
echo "  conda env remove -n scmina"
echo ""
echo "To update the environment:"
echo "  conda env update -f environment.yml"
echo ""
echo "Environment created. Activate with: conda activate scmina"

echo ""
echo "🎉 Environment setup complete!"
echo "Both Python and R are ready for scMINA development and Nextflow integration."
