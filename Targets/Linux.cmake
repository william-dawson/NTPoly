################################################################################
# Build file for a gcc, linux system.
set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_Fortran_COMPILER mpif90)
set(CMAKE_CXX_COMPILER mpicxx)

# Library Files
set(TOOLCHAIN_LIBS "-lblas")

# Release suggestions
set(CXX_TOOLCHAINFLAGS_RELEASE "-O3 -fopenmp")
set(F_TOOLCHAINFLAGS_RELEASE "-O3 -cpp -fopenmp")

# Debug suggestions
set(CXX_TOOLCHAINFLAGS_DEBUG "-O0 -fopenmp -Wall -Wextra \
    -pedantic -pedantic-errors -Wshadow -Wnull-dereference \
    -Wdouble-promotion -Woverloaded-virtual -Wmisleading-indentation \
    -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wuseless-cast \
    -Wsign-conversion -Wconversion -Wcast-align -Wold-style-cast")
set(F_TOOLCHAINFLAGS_DEBUG "-O0 -cpp -fcheck=all -Wall -Wextra -Werror \
    -pedantic -fimplicit-none -ffpe-trap=invalid,zero,overflow,underflow \
    -std=f2003")
