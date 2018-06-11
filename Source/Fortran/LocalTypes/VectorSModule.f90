!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling compressed vectors.
!! Compressed vectors are stored in two lists. The first is a list of indices,
!! the second a list of values.
MODULE VectorSModule
  USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: AddSparseVectors
  PUBLIC :: DotSparseVectors
  PUBLIC :: PairwiseMultiplyVectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE AddSparseVectors
     MODULE PROCEDURE AddSparseVectors_r
     MODULE PROCEDURE AddSparseVectors_c
  END INTERFACE
  INTERFACE DotSparseVectors
     MODULE PROCEDURE DotSparseVectors_r
     MODULE PROCEDURE DotSparseVectors_c
  END INTERFACE
  INTERFACE PairwiseMultiplyVectors
     MODULE PROCEDURE PairwiseMultiplyVectors_r
     MODULE PROCEDURE PairwiseMultiplyVectors_c
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define DATATYPE REAL(NTREAL)
#define AddSparseVectors AddSparseVectors_r
#define DotSparseVectors DotSparseVectors_r
#define PairwiseMultiplyVectors PairwiseMultiplyVectors_r
#include "includes/VectorSImpl.f90"
#undef AddSparseVectors
#undef DotSparseVectors
#undef PairwiseMultiplyVectors
#undef DATATYPE

#define DATATYPE COMPLEX(NTCOMPLEX)
#define AddSparseVectors AddSparseVectors_c
#define DotSparseVectors DotSparseVectors_c
#define PairwiseMultiplyVectors PairwiseMultiplyVectors_c
#include "includes/VectorSImpl.f90"
#undef AddSparseVectors
#undef DotSparseVectors
#undef PairwiseMultiplyVectors
#undef DATATYPE

END MODULE VectorSModule
