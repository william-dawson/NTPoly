################################################################################
# Target file for a Mac computer.
# Unfortunately, installing on the mac is a bit tricky. The problem is that
# when you install SWIG, it might link against the wrong version of python.
# That's why I've included custom set variables for the python dynamic library,
# which is determined by your choice of python executable.
set(CMAKE_SYSTEM_NAME Darwin)
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_Fortran_COMPILER mpif90)
set(CMAKE_CXX_COMPILER mpicxx)
set(PYTHON_EXECUTABLE python2)

execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c
"
from distutils.sysconfig import get_python_lib
path=get_python_lib(standard_lib=True)+\"/../../Python\"
print path"
    OUTPUT_VARIABLE PYTHON_LIBRARIES OUTPUT_STRIP_TRAILING_WHITESPACE)

set(CXX_TOOLCHAINFLAGS "-O3 -openmp -framework Accelerate -lgomp")
set(F_TOOLCHAINFLAGS "-O3 -cpp -fopenmp")
# Debug suggestions
#set(F_TOOLCHAINFLAGS "-fbounds-check -O0 -fexternal-blas -framework Accelerate -cpp -fopenmp -Wall -DPURE=")
