subdirs(Fortran)

if (NOT FORTRAN_ONLY)
  subdirs(C)
  subdirs(CPlusPlus)
  subdirs(Swig)
  subdirs(Wrapper)
endif()
