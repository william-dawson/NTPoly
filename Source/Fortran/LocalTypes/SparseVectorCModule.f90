!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling compressed vectors.
!! Compressed vectors are stored in two lists. The first is a list of indices,
!! the second a list of values.
!! This module can add two of those vectors together.
MODULE SparseVectorCModule
  USE DataTypesModule, ONLY : NTCOMPLEX, NTREAL

#define DATATYPE COMPLEX(NTCOMPLEX)
#include "includes/SparseVectorImplementation.f90"

END MODULE SparseVectorCModule
