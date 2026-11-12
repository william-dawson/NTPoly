# Genuine host-only NTPoly baseline: OpenACC directives are ignored and no
# GPU memory mode is enabled. Build the trs2_reference target only.
set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_Fortran_COMPILER mpif90)
set(CMAKE_CXX_COMPILER mpicxx)
# Export CPU fallback dependencies through the library target as well.
set(LAPACK_LIBRARIES lapack)
set(BLAS_LIBRARIES blas)
set(TOOLCHAIN_LIBS "-llapack -lblas")
set(F_TOOLCHAINFLAGS_RELEASE "-O3 -mp")
set(CXX_TOOLCHAINFLAGS_RELEASE "-O3")
set(F_TOOLCHAINFLAGS_DEBUG "-O0 -g -mp -Mbounds")
set(CXX_TOOLCHAINFLAGS_DEBUG "-O0 -g")
