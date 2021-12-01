################################################################################
# Target file for a Mac computer using anaconda. The important ingredient
# is modifying the find strategy.
set(CMAKE_SYSTEM_NAME Darwin)
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_Fortran_COMPILER mpif90)
set(CMAKE_CXX_COMPILER mpicxx)
set(PYTHON_EXECUTABLE python)
set(Python_FIND_STRATEGY LOCATION)

# Library Files
set(TOOLCHAIN_LIBS "-framework Accelerate -lgomp -L/usr/local/Cellar/scalapack/2.1.0_3/lib/ -lscalapack")

# Release Suggestions
set(CXX_TOOLCHAINFLAGS_RELEASE "-O3 -openmp")
set(F_TOOLCHAINFLAGS_RELEASE "-O3 -cpp -fopenmp -Wall")

# Debug suggestions
set(CXX_TOOLCHAINFLAGS_DEBUG "-O0")
set(F_TOOLCHAINFLAGS_DEBUG
  "-fcheck=all -O0 -fexternal-blas -cpp -Wall -DPURE=")
