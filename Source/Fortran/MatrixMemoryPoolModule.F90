!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling scratch memory for matrix multiplication.
!> The purpose of this module is to avoid having to allocate memory on the
!> heap during a matrix multiply, and to manage the underlying hash table.
MODULE MatrixMemoryPoolModule
  USE DataTypesModule, ONLY: NTREAL, NTCOMPLEX
  USE TripletModule, ONLY : Triplet_r, Triplet_c
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A memory pool datatype that can be reused for matrix matrix multiplication.
  !> this is to prevent excessive alloc/dealloc.
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
  !> this is to prevent excessive alloc/dealloc.
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
  PUBLIC :: ConstructMatrixMemoryPool
  PUBLIC :: DestructMatrixMemoryPool
  PUBLIC :: CheckMemoryPoolValidity
  PUBLIC :: SetPoolSparsity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE MatrixMemoryPool_lr
     MODULE PROCEDURE ConstructMatrixMemoryPool_lr
  END INTERFACE MatrixMemoryPool_lr
  INTERFACE MatrixMemoryPool_lc
     MODULE PROCEDURE ConstructMatrixMemoryPool_lc
  END INTERFACE MatrixMemoryPool_lc
  INTERFACE ConstructMatrixMemoryPool
     MODULE PROCEDURE ConstructMatrixMemoryPoolSub_lr
     MODULE PROCEDURE ConstructMatrixMemoryPoolSub_lc
  END INTERFACE ConstructMatrixMemoryPool
  INTERFACE DestructMatrixMemoryPool
     MODULE PROCEDURE DestructMatrixMemoryPool_lr
     MODULE PROCEDURE DestructMatrixMemoryPool_lc
  END INTERFACE DestructMatrixMemoryPool
  INTERFACE CheckMemoryPoolValidity
     MODULE PROCEDURE CheckMemoryPoolValidity_lr
     MODULE PROCEDURE CheckMemoryPoolValidity_lc
  END INTERFACE CheckMemoryPoolValidity
  INTERFACE SetPoolSparsity
     MODULE PROCEDURE SetPoolSparsity_lr
     MODULE PROCEDURE SetPoolSparsity_lc
  END INTERFACE SetPoolSparsity
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Subroutine wrapper for the constructor.
  SUBROUTINE ConstructMatrixMemoryPoolSub_lr(this, columns, rows, sparsity_in)
    !> The matrix to construct.
    TYPE(MatrixMemoryPool_lr), TARGET :: this
    !> Number of columns in the matrix.
    INTEGER, INTENT(IN) :: columns
    !> Number of rows in the matrix.
    INTEGER, INTENT(IN) :: rows
    !> Estimated sparsity (optional).
    REAL(NTREAL), INTENT(IN), OPTIONAL :: sparsity_in

#include "dense_includes/ConstructMatrixMemoryPool.f90"

  END SUBROUTINE ConstructMatrixMemoryPoolSub_lr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Subroutine wrapper for the constructor.
  SUBROUTINE ConstructMatrixMemoryPoolSub_lc(this, columns, rows, sparsity_in)
    !> The matrix to construct.
    TYPE(MatrixMemoryPool_lc), TARGET :: this
    !> Number of columns in the matrix.
    INTEGER, INTENT(IN) :: columns
    !> Number of rows in the matrix.
    INTEGER, INTENT(IN) :: rows
    !> Estimated sparsity (optional).
    REAL(NTREAL), INTENT(IN), OPTIONAL :: sparsity_in

#include "dense_includes/ConstructMatrixMemoryPool.f90"

  END SUBROUTINE ConstructMatrixMemoryPoolSub_lc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct Matrix Memory Pool object.
  FUNCTION ConstructMatrixMemoryPool_lr(columns, rows, sparsity_in) RESULT(this)
    !> The matrix to construct.
    TYPE(MatrixMemoryPool_lr), TARGET :: this
    !> Number of columns in the matrix.
    INTEGER, INTENT(IN) :: columns
    !> Number of rows in the matrix.
    INTEGER, INTENT(IN) :: rows
    !> Estimated sparsity (optional).
    REAL(NTREAL), INTENT(IN), OPTIONAL :: sparsity_in

#include "dense_includes/ConstructMatrixMemoryPool.f90"

  END FUNCTION ConstructMatrixMemoryPool_lr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct Matrix Memory Pool object.
  FUNCTION ConstructMatrixMemoryPool_lc(columns, rows, sparsity_in) RESULT(this)
    !> The matrix to construct.
    TYPE(MatrixMemoryPool_lc), TARGET :: this
    !> Number of columns in the matrix.
    INTEGER, INTENT(IN) :: columns
    !> Number of rows in the matrix.
    INTEGER, INTENT(IN) :: rows
    !> Estimated sparsity (optional).
    REAL(NTREAL), INTENT(IN), OPTIONAL :: sparsity_in

#include "dense_includes/ConstructMatrixMemoryPool.f90"

  END FUNCTION ConstructMatrixMemoryPool_lc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A destructor for a matrix memory pool
  PURE SUBROUTINE DestructMatrixMemoryPool_lr(this)
    !> The matrix being destructed.
    TYPE(MatrixMemoryPool_lr), INTENT(INOUT) :: this

#include "dense_includes/DestructMatrixMemoryPool.f90"

  END SUBROUTINE DestructMatrixMemoryPool_lr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A destructor for a matrix memory pool
  PURE SUBROUTINE DestructMatrixMemoryPool_lc(this)
    !> The matrix being destructed.
    TYPE(MatrixMemoryPool_lc), INTENT(INOUT) :: this

#include "dense_includes/DestructMatrixMemoryPool.f90"

  END SUBROUTINE DestructMatrixMemoryPool_lc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Checks if a given memory pool has been validly allocated to handle
  !> the given parameters.
  PURE FUNCTION CheckMemoryPoolValidity_lr(this, columns, rows) RESULT(isvalid)
    !> The memory pool to check.
    TYPE(MatrixMemoryPool_lr), INTENT(in) :: this
    !> Number of columns in the matrix.
    INTEGER, INTENT(IN) :: columns
    !> Number of rows in the matrix.
    INTEGER, INTENT(IN) :: rows
    !> true if the memory pool is valid.
    LOGICAL :: isvalid

#include "dense_includes/CheckMemoryPoolValidity.f90"

  END FUNCTION CheckMemoryPoolValidity_lr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Checks if a given memory pool has been validly allocated to handle
  !> Checks if a given memory pool has been validly allocated to handle
  !> the given parameters.
  PURE FUNCTION CheckMemoryPoolValidity_lc(this, columns, rows) RESULT(isvalid)
    !> The memory pool to check.
    TYPE(MatrixMemoryPool_lc), INTENT(in) :: this
    !> Number of columns in the matrix.
    INTEGER, INTENT(IN) :: columns
    !> Number of rows in the matrix.
    INTEGER, INTENT(IN) :: rows
    !> true if the memory pool is valid.
    LOGICAL :: isvalid

#include "dense_includes/CheckMemoryPoolValidity.f90"

  END FUNCTION CheckMemoryPoolValidity_lc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets the expected sparsity of the matrix, which helps with hashing.
  SUBROUTINE SetPoolSparsity_lr(this,sparsity)
    !> The memory pool to set the sparsity of.
    TYPE(MatrixMemoryPool_lr), INTENT(INOUT), TARGET :: this
    !> The sparsity value.
    REAL(NTREAL), INTENT(IN) :: sparsity

#include "dense_includes/SetPoolSparsity.f90"

  END SUBROUTINE SetPoolSparsity_lr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets the expected sparsity of the matrix, which helps with hashing.
  SUBROUTINE SetPoolSparsity_lc(this,sparsity)
    !> The memory pool to set the sparsity of.
    TYPE(MatrixMemoryPool_lc), INTENT(INOUT), TARGET :: this
    !> The sparsity value.
    REAL(NTREAL), INTENT(IN) :: sparsity

#include "dense_includes/SetPoolSparsity.f90"

  END SUBROUTINE SetPoolSparsity_lc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE MatrixMemoryPoolModule
