################################################################################
# Build file for a gcc, linux system.
set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_Fortran_COMPILER mpif90)
set(CMAKE_CXX_COMPILER mpicxx)

# Release suggestions
set(CXX_TOOLCHAINFLAGS_RELEASE "-O3 -openmp -lgomp -fPIC -llapack")
set(F_TOOLCHAINFLAGS_RELEASE "-O3 -cpp -openmp -fPIC")

# Debug suggestions
set(CXX_TOOLCHAINFLAGS_DEBUG "-O0 -openmp -fPIC -llapack")
set(F_TOOLCHAINFLAGS_DEBUG "-O0 -cpp -fPIC")
