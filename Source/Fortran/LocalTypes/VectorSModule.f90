!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling compressed vectors.
!! Compressed vectors are stored in two lists. The first is a list of indices,
!! the second a list of values.
!! This module can add two of those vectors together.
MODULE VectorSRModule
  USE DataTypesModule, ONLY : NTREAL

#define DATATYPE REAL(NTREAL)

#include "includes/VectorSImpl.f90"

#undef DATATYPE

END MODULE VectorSRModule

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling compressed vectors.
!! Compressed vectors are stored in two lists. The first is a list of indices,
!! the second a list of values.
!! This module can add two of those vectors together.
MODULE VectorSCModule
  USE DataTypesModule, ONLY : NTCOMPLEX, NTREAL

#define DATATYPE COMPLEX(NTCOMPLEX)

#include "includes/VectorSImpl.f90"

#undef DATATYPE

END MODULE VectorSCModule
