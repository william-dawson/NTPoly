################################################################################
# Build file for an intel compiler system.
set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_C_COMPILER mpiicc)
set(CMAKE_Fortran_COMPILER mpiifort)
set(CMAKE_CXX_COMPILER mpiicpc)

set(CXX_TOOLCHAINFLAGS "-O3 -qopenmp -lgomp -fPIC")
set(F_TOOLCHAINFLAGS "-O3 -fpp -qopenmp -fPIC")
# Debug suggestions
#set(CXX_TOOLCHAINFLAGS "-O0 -qopenmp -lgomp -fPIC")
#set(F_TOOLCHAINFLAGS "-check bounds -O0 -fpp -qopenmp -fPIC -DPURE=")
