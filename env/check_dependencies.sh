#!/bin/bash
# Dependency checker and installer for LAO-STO project

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

CONDA_ENV="SC"
MISSING_DEPS=()

echo "========================================="
echo "LAO-STO Dependency Checker"
echo "========================================="

# Check Intel Fortran Compiler (ifx)
echo -e "\n[1/3] Checking Intel Fortran Compiler (ifx)..."
if command -v ifx &> /dev/null; then
    VERSION=$(ifx --version | head -n1)
    echo -e "${GREEN}✓${NC} ifx found: $VERSION"
else
    echo -e "${RED}✗${NC} ifx not found"
    MISSING_DEPS+=("ifx")
fi

# Check pfUnit
echo -e "\n[2/3] Checking pfUnit..."
if [ -n "$PFUNIT_DIR" ] && [ -d "$PFUNIT_DIR" ]; then
    echo -e "${GREEN}✓${NC} pfUnit found at: $PFUNIT_DIR"
else
    echo -e "${RED}✗${NC} pfUnit not found (set \$PFUNIT_DIR or install)"
    MISSING_DEPS+=("pfunit")
fi

# Check Conda environment
echo -e "\n[3/3] Checking Conda environment '$CONDA_ENV'..."
if command -v conda &> /dev/null; then
    if conda env list | grep -q "^$CONDA_ENV "; then
        echo -e "${GREEN}✓${NC} Conda environment '$CONDA_ENV' exists"
        
        # Activate and check Python packages
        eval "$(conda shell.bash hook)"
        conda activate $CONDA_ENV 2>/dev/null || true
        
        echo -e "\nChecking Python packages in '$CONDA_ENV'..."
        PYTHON_DEPS=("numpy" "matplotlib" "scipy" "sympy" "f90nml" "fortranformat" "colorcet" "sklearn" "yaml")
        
        for pkg in "${PYTHON_DEPS[@]}"; do
            if python -c "import ${pkg}" 2>/dev/null; then
                echo -e "  ${GREEN}✓${NC} $pkg"
            else
                echo -e "  ${RED}✗${NC} $pkg"
                MISSING_DEPS+=("python:$pkg")
            fi
        done
    else
        echo -e "${RED}✗${NC} Conda environment '$CONDA_ENV' not found"
        MISSING_DEPS+=("conda_env")
    fi
else
    echo -e "${RED}✗${NC} Conda not found"
    MISSING_DEPS+=("conda")
fi

# Summary
echo -e "\n========================================="
if [ ${#MISSING_DEPS[@]} -eq 0 ]; then
    echo -e "${GREEN}All dependencies satisfied!${NC}"
    exit 0
else
    echo -e "${YELLOW}Missing dependencies:${NC}"
    for dep in "${MISSING_DEPS[@]}"; do
        echo "  - $dep"
    done
    
    echo -e "\n${YELLOW}Installation suggestions:${NC}"
    
    for dep in "${MISSING_DEPS[@]}"; do
        case $dep in
            ifx)
                echo -e "\n• Intel Fortran Compiler (ifx):"
                echo "  Download from: https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html"
                echo "  Or install Intel oneAPI Base & HPC Toolkit (recommended)."
                ;;
            pfunit)
                echo -e "\n• pfUnit:"
                echo "  git clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git"
                echo "  cd pFUnit && mkdir build && cd build"
                echo "  cmake .. && make && make install"
                echo "  export PFUNIT_DIR=<path_to_pfunit_install>"
                ;;
            conda_env)
                echo -e "\n• Conda environment '$CONDA_ENV':"
                echo "  conda create -n $CONDA_ENV python=3.10"
                echo "  conda activate $CONDA_ENV"
                echo "  pip install numpy matplotlib scipy sympy f90nml fortranformat colorcet scikit-learn pyyaml"
                ;;
            python:*)
                PKG="${dep#python:}"
                echo "  pip install $PKG"
                ;;
        esac
    done
    
    exit 1
fi
