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
  !> Construct Matrix Memory Pool object.
  !> @param[out] this a constructed Matrix Memory Pool object.
  !> @param[in] columns number of columns in the matrix.
  !> @param[in] rows number of rows in the matrix.
  !> @param[in] sparsity_in estimated sparsity (optional).
  FUNCTION ConstructMatrixMemoryPool_lr(columns, rows, sparsity_in) RESULT(this)
    !! Parameters
    TYPE(MatrixMemoryPool_lr), TARGET :: this
    INTEGER(kind=c_int), INTENT(IN) :: columns
    INTEGER(kind=c_int), INTENT(IN) :: rows
    REAL(NTREAL), INTENT(IN), OPTIONAL :: sparsity_in

    INCLUDE "includes/ConstructMatrixMemoryPool.f90"

  END FUNCTION ConstructMatrixMemoryPool_lr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A destructor for a matrix memory pool
  !> @param[inout] this the matrix being destructed.
  PURE SUBROUTINE DestructMatrixMemoryPool_lr(this)
    !! Parameters
    TYPE(MatrixMemoryPool_lr), INTENT(INOUT) :: this

    INCLUDE "includes/DestructMatrixMemoryPool.f90"

  END SUBROUTINE DestructMatrixMemoryPool_lr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Checks if a given memory pool has been validly allocated to handle
  !! the given parameters.
  !> @param[in] this the memory pool to check.
  !> @param[in] columns number of columns in the matrix.
  !> @param[in] rows number of rows in the matrix.
  !> @return true if the memory pool is valid.
  PURE FUNCTION CheckMemoryPoolValidity_lr(this, columns, rows) RESULT(isvalid)
    !! Parameters
    TYPE(MatrixMemoryPool_lr), INTENT(in) :: this
    INTEGER, INTENT(IN) :: columns
    INTEGER, INTENT(IN) :: rows
    LOGICAL :: isvalid

    INCLUDE "includes/CheckMemoryPoolValidity.f90"

  END FUNCTION CheckMemoryPoolValidity_lr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets the expected sparsity of the matrix, which helps with hashing.
  !! @param[inout] this the memory pool to set the sparsity of.
  !! @param[in] sparsity the sparsity value.
  SUBROUTINE SetPoolSparsity_lr(this,sparsity)
    !! Parameters
    TYPE(MatrixMemoryPool_lr), INTENT(INOUT), TARGET :: this
    REAL(NTREAL), INTENT(IN) :: sparsity

    INCLUDE "includes/SetPoolSparsity.f90"

  END SUBROUTINE SetPoolSparsity_lr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct Matrix Memory Pool object.
  !> @param[out] this a constructed Matrix Memory Pool object.
  !> @param[in] columns number of columns in the matrix.
  !> @param[in] rows number of rows in the matrix.
  !> @param[in] sparsity_in estimated sparsity (optional).
  FUNCTION ConstructMatrixMemoryPool_lc(columns, rows, sparsity_in) RESULT(this)
    !! Parameters
    TYPE(MatrixMemoryPool_lc), TARGET :: this
    INTEGER(kind=c_int), INTENT(IN) :: columns
    INTEGER(kind=c_int), INTENT(IN) :: rows
    REAL(NTREAL), INTENT(IN), OPTIONAL :: sparsity_in

    INCLUDE "includes/ConstructMatrixMemoryPool.f90"

  END FUNCTION ConstructMatrixMemoryPool_lc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A destructor for a matrix memory pool
  !> @param[inout] this the matrix being destructed.
  PURE SUBROUTINE DestructMatrixMemoryPool_lc(this)
    !! Parameters
    TYPE(MatrixMemoryPool_lc), INTENT(INOUT) :: this

    INCLUDE "includes/DestructMatrixMemoryPool.f90"

  END SUBROUTINE DestructMatrixMemoryPool_lc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Checks if a given memory pool has been validly allocated to handle
  !! the given parameters.
  !> @param[in] this the memory pool to check.
  !> @param[in] columns number of columns in the matrix.
  !> @param[in] rows number of rows in the matrix.
  !> @return true if the memory pool is valid.
  PURE FUNCTION CheckMemoryPoolValidity_lc(this, columns, rows) RESULT(isvalid)
    !! Parameters
    TYPE(MatrixMemoryPool_lc), INTENT(in) :: this
    INTEGER, INTENT(IN) :: columns
    INTEGER, INTENT(IN) :: rows
    LOGICAL :: isvalid

    INCLUDE "includes/CheckMemoryPoolValidity.f90"

  END FUNCTION CheckMemoryPoolValidity_lc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets the expected sparsity of the matrix, which helps with hashing.
  !! @param[inout] this the memory pool to set the sparsity of.
  !! @param[in] sparsity the sparsity value.
  SUBROUTINE SetPoolSparsity_lc(this,sparsity)
    !! Parameters
    TYPE(MatrixMemoryPool_lc), INTENT(INOUT), TARGET :: this
    REAL(NTREAL), INTENT(IN) :: sparsity

    INCLUDE "includes/SetPoolSparsity.f90"

  END SUBROUTINE SetPoolSparsity_lc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE MatrixMemoryPoolModule
