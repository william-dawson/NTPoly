################################################################################
# Build file for an intel mkl system.
set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_C_COMPILER mpiicc)
set(CMAKE_Fortran_COMPILER mpiifort)
set(CMAKE_CXX_COMPILER mpiicpc)

set(PYTHON_INCLUDE_PATH "/usr/include/python2.6/")
#set(CXX_TOOLCHAINFLAGS "-qopenmp -lgomp -fPIC")
#set(F_TOOLCHAINFLAGS "-check bounds -O0 -fpp -qopenmp -fPIC")
set(CXX_TOOLCHAINFLAGS "-O3 -qopenmp -lgomp -fPIC")
set(F_TOOLCHAINFLAGS "-O3 -fpp -qopenmp -fPIC")
