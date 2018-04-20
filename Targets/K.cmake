################################################################################
# Build file for the K Computer.
set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_C_COMPILER mpifccpx)
set(CMAKE_Fortran_COMPILER mpifrtpx)
set(CMAKE_CXX_COMPILER mpiFCCpx)

# Release Suggestions
set(CXX_TOOLCHAINFLAGS_RELEASE
    "-Kfast,-Kparallel,openmp,optmsg=2 --linkfortran")
set(F_TOOLCHAINFLAGS_RELEASE "-Kfast,-Kparallel,openmp,optmsg=2 -Cpp")

# Debug Suggestions
set(CXX_TOOLCHAINFLAGS_DEBUG "--linkfortran")
set(F_TOOLCHAINFLAGS_DEBUG "-Cpp")
