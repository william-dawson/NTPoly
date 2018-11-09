################################################################################
# Target file for a Mac computer.
# Unfortunately, installing on the mac is a bit tricky. The problem is that
# when you install SWIG, it might link against the wrong version of python.
# That's why I've included custom set variables for the python executable.
set(CMAKE_SYSTEM_NAME Darwin)
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_Fortran_COMPILER mpif90)
set(CMAKE_CXX_COMPILER mpicxx)
set(PYTHON_EXECUTABLE python3)

# Library Files
set(TOOLCHAIN_LIBS "-framework Accelerate -lgomp")

# Release Suggestions
set(CXX_TOOLCHAINFLAGS_RELEASE "-O3 -openmp")
set(F_TOOLCHAINFLAGS_RELEASE "-O3 -cpp -fopenmp")

# Debug suggestions
set(CXX_TOOLCHAINFLAGS_DEBUG "-O0")
set(F_TOOLCHAINFLAGS_DEBUG
  "-fbounds-check -O0 -fexternal-blas -cpp -Wall -DPURE=")
