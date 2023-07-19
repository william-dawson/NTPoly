set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_C_COMPILER mpiicc)
set(CMAKE_Fortran_COMPILER mpiifx)
set(CMAKE_CXX_COMPILER mpiicpc)

# Library Files
set(TOOLCHAIN_LIBS "-qmkl=sequential")

# Release suggestions
set(CXX_TOOLCHAINFLAGS_RELEASE "-O3 -qopenmp -lgomp")
set(F_TOOLCHAINFLAGS_RELEASE "-O3 -fpp -qopenmp")

# Debug suggestions
set(CXX_TOOLCHAINFLAGS_DEBUG "-O0 -qopenmp -lgomp")
set(F_TOOLCHAINFLAGS_DEBUG "-check bounds -O0 -g -fpp -traceback "
                           "-qopenmp -DPURE=")
