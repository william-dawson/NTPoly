!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling scratch memory for matrix multiplication.
!! The purpose of this module is to avoid having to allocate memory on the
!! heap during a matrix multiply, and to manage the underlying hash table.
MODULE MatrixMemoryPoolRModule
  USE DataTypesModule, ONLY: NTREAL
  USE TripletRModule, ONLY : Triplet_r

#define DATATYPE REAL(NTREAL)
#define TTYPE Triplet_r
#define MPOOLTYPE MatrixMemoryPool_r

#include "includes/MatrixMemoryPoolImpl.f90"

#undef MPOOLTYPE
#undef TTYPE
#undef DATATYPE

END MODULE MatrixMemoryPoolRModule

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

#include "includes/MatrixMemoryPoolImpl.f90"

#undef MPOOLTYPE
#undef TTYPE
#undef DATATYPE

END MODULE MatrixMemoryPoolCModule
