!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling scratch memory for matrix multiplication.
!! The purpose of this module is to avoid having to allocate memory on the
!! heap during a matrix multiply, and to manage the underlying hash table.
MODULE MatrixMemoryPoolModule
  USE DataTypesModule, ONLY: NTREAL
  USE TripletModule, ONLY : Triplet_t

#define DATATYPE REAL(NTREAL)
#define TTYPE Triplet_t
#define MPOOLTYPE MatrixMemoryPool_t

#include "includes/MatrixMemoryPoolImplementation.f90"

#undef MPOOLTYPE
#undef TTYPE
#undef DATATYPE

END MODULE MatrixMemoryPoolModule

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling scratch memory for matrix multiplication.
!! The purpose of this module is to avoid having to allocate memory on the
!! heap during a matrix multiply, and to manage the underlying hash table.
MODULE MatrixMemoryPoolCModule
  USE DataTypesModule, ONLY: NTCOMPLEX, NTREAL
  USE TripletCModule, ONLY : Triplet_C

#define DATATYPE COMPLEX(NTCOMPLEX)
#define TTYPE Triplet_c
#define MPOOLTYPE MatrixMemoryPool_c

#include "includes/MatrixMemoryPoolImplementation.f90"

#undef MPOOLTYPE
#undef TTYPE
#undef DATATYPE

END MODULE MatrixMemoryPoolCModule
