################################################################################
## Basic Code
subdirs(Fortran)

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
