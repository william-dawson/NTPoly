!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling compressed vectors.
!! Compressed vectors are stored in two lists. The first is a list of indices,
!! the second a list of values.
!! This module can add two of those vectors together.
MODULE SparseVectorModule
  USE DataTypesModule, ONLY : NTREAL

#define DATATYPE REAL(NTREAL)
#include "includes/SparseVectorImplementation.f90"

END MODULE SparseVectorModule
