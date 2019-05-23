################################################################################
# Build file for the archer system.
# Using the default cray compilers
set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_C_COMPILER cc) 
set(CMAKE_Fortran_COMPILER ftn)
set(CMAKE_CXX_COMPILER CC)

# Library Files
#set(TOOLCHAIN_LIBS "-lblas")

# Release suggestions
set(CXX_TOOLCHAINFLAGS_RELEASE "-O3 -dynamic")
set(F_TOOLCHAINFLAGS_RELEASE "-O3 -dynamic")

# Debug suggestions
set(CXX_TOOLCHAINFLAGS_DEBUG "-dynamic")
set(F_TOOLCHAINFLAGS_DEBUG "-dynamic -DPURE=")

set(TARGET_SUPPORTS_SHARED_LIBS Yes)
