# About f2pyfoil v1.0
f2pyfoil is a project aimed at creating a Python interface for the XFOIL airfoil analysis software using Fortran. The goal of this project is to provide a high-performance and simple interface for XFOIL that can be easily integrated into Python workflows.

Currently the only part of XFOIL that is implemented in f2pyfoil is the panel method solver, which is used to calculate the flow around an airfoil. The panel method solver is a key component of XFOIL and is responsible for calculating the pressure distribution and aerodynamic coefficients for a given airfoil geometry and flow conditions. Future plans aim to implement XFOIL's remaining features to provide a complete version of XFOIL.

This repository will be exclusively for the development of the Fortran side of the code and will provide no Python code. The Fortran code will be compiled into a shared library that can be imported and used in Python using the f2py tool. The Python wrapper code will be developed in a separate repository that will depend on this one.

# Future Plans (v 1.1)
- [x] Add boundary layer data retrieval functions
- [x] Add parameterised airfoil geometry definitions (PARSEC, CST, etc.)
- [x] Add Laitone and Prandtl-Glauert compressibility corrections

- [ ] Add support for other analysis routines (Surface-speed (QDES), Inverse (MDES) and Geometry (GDES))


# f2pyfoil build guide

The build has been tested on a Windows 11 64-bit system with Pytohn 3.11 and gfortran 15.2.0. Static linking for some gfortran libraries is required to avoid runtime errors. The build process uses the `uv` package manager to create an isolated environment and manage dependencies. Conda has known issues for building Fortran code on Windows particularly with gfortran.

## Requirements

| Requirement | Version | Notes |
|-------------|---------|-------|
| Python | 3.9 or later | 3.11 recommended |
| gfortran | 9.0 or later | Must be on system PATH |
| uv | latest | Package manager |

## Help with installing requirements
1. **Install gfortran**: NumPy provides an installation guide for gfortran on Windows: https://numpy.org/doc/2.1/f2py/windows/msys2.html. If gfortran is successfully installed and visible on the PATH running the command "gfortran --version" should print the version number in the terminal.

2. **Install uv**: The `uv` package manager can be installed following the guide on its official website: https://docs.astral.sh/uv/getting-started/installation/. Winget was used for the test build

3. **Python Environment**: It is recommended to use a virtual environment to manage dependencies. You can create and activate a virtual environment using `uv` following its guide. For example in PowerShell you can run the following command to setup the virtual environment with Python 3.11:
   ```
    uv venv --python 3.11
    .venv\Scripts\activate
    uv pip install meson ninja meson-python numpy
   ```
## Build Instructions
1. **Clone the Repository**: Clone the f2pyfoil repository to your local machine.
   ```
   git clone https://github.com/C-Williams-Projects/f2pyfoil
    cd f2pyfoil
    ```
    [activate your environment here, if not already activated]
    ```
    uv pip install --no-buld-isolation -e .
    ```

2. **Quick Check**: After the installation, you can run a quick check to ensure that the package is installed correctly. This should print a long list of information abour the xfoil module, including the available functions and classes.
   ```
    uv run python -c "import xfoil; print(dir(xfoil))"
   ```

## Full Test
To run the full test suite, you can use the test.py file included in the repository. This will run a series of tests to verify that the wrapper is functioning correctly.

The reference data is taken from a single precision version of XFOIL. Additionally, it is a precompiled version of XFOIL for Windows so discrepancies in the data are expected when running on other platforms or with different versions of gfortran. The test suite is designed to check for consistency with the reference data, so some differences in the results may be acceptable as long as they are within a reasonable range.
