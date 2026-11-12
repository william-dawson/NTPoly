################################################################################
## Basic Code
add_subdirectory(Fortran)
if(NTPOLY_ENABLE_CUSOLVERMP)
  add_subdirectory(GPU)
endif()

## C++ Bindings
if (NOT FORTRAN_ONLY)
  subdirs(C)
  subdirs(CPlusPlus)
  subdirs(Wrapper)
endif()

## Python bindings built with swig
if (NOT NOSWIG)
  subdirs(Swig)
else()
  message(WARNING "Swig is disabled. No python bindings with be generated.")
endif()
