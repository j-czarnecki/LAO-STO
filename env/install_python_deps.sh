#!/bin/bash
# Auto-installer for LAO-STO Python dependencies

set -e

CONDA_ENV="SC"
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo "========================================="
echo "LAO-STO Python Dependency Installer"
echo "========================================="

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo -e "${RED}Error: Conda not found. Please install Miniconda or Anaconda first.${NC}"
    exit 1
fi

# Initialize conda for bash
eval "$(conda shell.bash hook)"

# Check if environment exists, create if not
if ! conda env list | grep -q "^$CONDA_ENV "; then
    echo -e "${YELLOW}Creating conda environment '$CONDA_ENV'...${NC}"
    conda create -n $CONDA_ENV python=3.10 -y
fi

# Activate environment
echo -e "\n${GREEN}Activating conda environment '$CONDA_ENV'...${NC}"
conda activate $CONDA_ENV

# Install Python packages
echo -e "\n${GREEN}Installing Python packages...${NC}"
PYTHON_DEPS=(
    "numpy"
    "matplotlib"
    "scipy"
    "sympy"
    "f90nml"
    "fortranformat"
    "colorcet"
    "scikit-learn"
    "pyyaml"
    "ipython"
    "jupyter"
    "pandas"
)

for pkg in "${PYTHON_DEPS[@]}"; do
    echo -e "\nInstalling $pkg..."
    pip install "$pkg" --upgrade
done

echo -e "\n${GREEN}=========================================${NC}"
echo -e "${GREEN}Python dependencies installed successfully!${NC}"
echo -e "${GREEN}=========================================${NC}"
echo -e "\nTo activate the environment, run:"
echo -e "  ${YELLOW}conda activate $CONDA_ENV${NC}"
echo -e "\nNote: Intel Fortran Compiler (ifx) and pfUnit must be installed manually."
