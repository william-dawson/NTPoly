################################################################################
# Build file for the K Computer.
set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_C_COMPILER mpifccpx)
set(CMAKE_Fortran_COMPILER mpifrtpx)
set(CMAKE_CXX_COMPILER mpiFCCpx)
SET(CMAKE_Fortran_MODDIR_FLAG "-M")

set(TOOLCHAIN_LIBS "")

# Release Suggestions
set(CXX_TOOLCHAINFLAGS_RELEASE
    "-Kfast,parallel,openmp --linkfortran")
set(F_TOOLCHAINFLAGS_RELEASE "-Kfast,parallel,openmp -Cpp")

# Debug Suggestions
set(CXX_TOOLCHAINFLAGS_DEBUG "--linkfortran")
set(F_TOOLCHAINFLAGS_DEBUG "-Cpp")
