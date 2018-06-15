!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling scratch memory for matrix multiplication.
!! The purpose of this module is to avoid having to allocate memory on the
!! heap during a matrix multiply, and to manage the underlying hash table.
MODULE MatrixMemoryPoolModule
  USE DataTypesModule, ONLY: NTREAL, NTCOMPLEX
  USE TripletModule, ONLY : Triplet_r, Triplet_c
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A memory pool datatype that can be reused for matrix matrix multiplication.
  !! this is to prevent excessive alloc/dealloc.
  TYPE, PUBLIC :: MatrixMemoryPool_lr
     PRIVATE
     !> Shape of matrix: columns
     INTEGER, PUBLIC :: columns
     !> Shape of matrix: rows
     INTEGER, PUBLIC :: rows
     !> storage for actual values added to the matrix.
     TYPE(Triplet_r), DIMENSION(:), ALLOCATABLE, PUBLIC :: pruned_list
     !> storage for potential values added to the matrix.
     REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: value_array
     !> true if an element has been pushed to this part of the matrix.
     LOGICAL, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: dirty_array
     !> Storage space for indices, hashed.
     INTEGER, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: hash_index
     !> Internal storage space for amount of items added to a bucket.
     INTEGER, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: inserted_per_bucket
     !> Size of the buckets.
     INTEGER, PUBLIC :: hash_size
  END TYPE MatrixMemoryPool_lr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A memory pool datatype that can be reused for matrix matrix multiplication.
  !! this is to prevent excessive alloc/dealloc.
  TYPE, PUBLIC :: MatrixMemoryPool_lc
     PRIVATE
     !> Shape of matrix: columns
     INTEGER, PUBLIC :: columns
     !> Shape of matrix: rows
     INTEGER, PUBLIC :: rows
     !> storage for actual values added to the matrix.
     TYPE(Triplet_c), DIMENSION(:), ALLOCATABLE, PUBLIC :: pruned_list
     !> storage for potential values added to the matrix.
     COMPLEX(NTCOMPLEX), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: value_array
     !> true if an element has been pushed to this part of the matrix.
     LOGICAL, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: dirty_array
     !> Storage space for indices, hashed.
     INTEGER, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: hash_index
     !> Internal storage space for amount of items added to a bucket.
     INTEGER, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: inserted_per_bucket
     !> Size of the buckets.
     INTEGER, PUBLIC :: hash_size
  END TYPE MatrixMemoryPool_lc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: DestructMatrixMemoryPool
  PUBLIC :: CheckMemoryPoolValidity
  PUBLIC :: SetPoolSparsity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE MatrixMemoryPool_lr
     MODULE PROCEDURE ConstructMatrixMemoryPool_lr
  END INTERFACE
  INTERFACE MatrixMemoryPool_lc
     MODULE PROCEDURE ConstructMatrixMemoryPool_lc
  END INTERFACE
  INTERFACE DestructMatrixMemoryPool
     MODULE PROCEDURE DestructMatrixMemoryPool_lr
     MODULE PROCEDURE DestructMatrixMemoryPool_lc
  END INTERFACE
  INTERFACE CheckMemoryPoolValidity
     MODULE PROCEDURE CheckMemoryPoolValidity_lr
     MODULE PROCEDURE CheckMemoryPoolValidity_lc
  END INTERFACE
  INTERFACE SetPoolSparsity
     MODULE PROCEDURE SetPoolSparsity_lr
     MODULE PROCEDURE SetPoolSparsity_lc
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define DATATYPE REAL(NTREAL)
#define MPOOLTYPE MatrixMemoryPool_lr
#define ConstructMatrixMemoryPool ConstructMatrixMemoryPool_lr
#define DestructMatrixMemoryPool DestructMatrixMemoryPool_lr
#define CheckMemoryPoolValidity CheckMemoryPoolValidity_lr
#define SetPoolSparsity SetPoolSparsity_lr

#include "includes/MatrixMemoryPoolImpl.f90"

#undef ConstructMatrixMemoryPool
#undef DestructMatrixMemoryPool
#undef CheckMemoryPoolValidity
#undef SetPoolSparsity
#undef MPOOLTYPE
#undef TTYPE
#undef DATATYPE

#define DATATYPE COMPLEX(NTCOMPLEX)
#define MPOOLTYPE MatrixMemoryPool_lc
#define ConstructMatrixMemoryPool ConstructMatrixMemoryPool_lc
#define DestructMatrixMemoryPool DestructMatrixMemoryPool_lc
#define CheckMemoryPoolValidity CheckMemoryPoolValidity_lc
#define SetPoolSparsity SetPoolSparsity_lc

#include "includes/MatrixMemoryPoolImpl.f90"
#undef ConstructMatrixMemoryPool
#undef DestructMatrixMemoryPool
#undef CheckMemoryPoolValidity
#undef SetPoolSparsity
#undef MPOOLTYPE
#undef TTYPE
#undef DATATYPE

END MODULE MatrixMemoryPoolModule
