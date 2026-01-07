# Dependency Management Scripts

This directory contains scripts to check and install dependencies for the LAO-STO project.

## Scripts

### 1. `check_dependencies.sh`
Checks if all required dependencies are installed:
- Intel Fortran Compiler (ifx)
- pfUnit (for unit testing)
- Python environment (conda environment named "SC")
- Python packages: numpy, matplotlib, scipy, sympy, f90nml, fortranformat, colorcet, scikit-learn, pyyaml

**Usage:**
```bash
bash check_dependencies.sh
```

**Exit codes:**
- 0: All dependencies satisfied
- 1: Missing dependencies (with installation suggestions)

### 2. `install_python_deps.sh`
Automatically installs Python dependencies in the conda environment "SC".
Creates the environment if it doesn't exist.

**Usage:**
```bash
bash install_python_deps.sh
```

**Note:** This script only installs Python packages. Intel Fortran Compiler and pfUnit must be installed manually.
