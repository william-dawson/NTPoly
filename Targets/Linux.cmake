################################################################################
# Build file for a gcc, linux system.
set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_Fortran_COMPILER mpif90)
set(CMAKE_CXX_COMPILER mpicxx)

set(CXX_TOOLCHAINFLAGS "-O3 -openmp -lgomp -fPIC")
set(F_TOOLCHAINFLAGS "-O3 -cpp -openmp -fPIC")
