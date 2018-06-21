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
  TYPE, ABSTRACT, PUBLIC :: MatrixMemoryPool_l
     !> Shape of matrix: columns
     INTEGER, PUBLIC :: columns
     !> Shape of matrix: rows
     INTEGER, PUBLIC :: rows
     !> true if an element has been pushed to this part of the matrix.
     LOGICAL, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: dirty_array
     !> Storage space for indices, hashed.
     INTEGER, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: hash_index
     !> Internal storage space for amount of items added to a bucket.
     INTEGER, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: inserted_per_bucket
     !> Size of the buckets.
     INTEGER, PUBLIC :: hash_size
  CONTAINS
     PROCEDURE(Construct_base), DEFERRED :: Init
     PROCEDURE(Destruct_base), DEFERRED :: Destruct
     PROCEDURE(CheckValidity_base), DEFERRED :: CheckValidity
     PROCEDURE(SetSparsity_base), DEFERRED :: SetSparsity
  END TYPE MatrixMemoryPool_l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE, EXTENDS(MatrixMemoryPool_l), PUBLIC :: MatrixMemoryPool_lr
     !> storage for actual values added to the matrix.
     TYPE(Triplet_r), DIMENSION(:), ALLOCATABLE, PUBLIC :: pruned_list
     !> storage for potential values added to the matrix.
     REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: value_array
  CONTAINS
     PROCEDURE :: Init => Construct_r
     PROCEDURE :: Destruct => Destruct_r
     PROCEDURE :: CheckValidity => CheckValidity_r
     PROCEDURE :: SetSparsity => SetSparsity_r
  END TYPE MatrixMemoryPool_lr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A memory pool datatype that can be reused for matrix matrix multiplication.
  !! this is to prevent excessive alloc/dealloc.
  TYPE, EXTENDS(MatrixMemoryPool_l), PUBLIC :: MatrixMemoryPool_lc
     TYPE(Triplet_c), DIMENSION(:), ALLOCATABLE, PUBLIC :: pruned_list
     !> storage for potential values added to the matrix.
     COMPLEX(NTCOMPLEX), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: value_array
  CONTAINS
     PROCEDURE :: Init => Construct_c
     PROCEDURE :: Destruct => Destruct_c
     PROCEDURE :: CheckValidity => CheckValidity_c
     PROCEDURE :: SetSparsity => SetSparsity_c
  END TYPE MatrixMemoryPool_lc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ABSTRACT INTERFACE
    SUBROUTINE Construct_base(this, columns, rows, sparsity_in)
      USE DataTypesModule, ONLY : NTREAL
      USE ISO_C_BINDING, ONLY : c_int
      IMPORT :: MatrixMemoryPool_l
      IMPLICIT NONE
      CLASS(MatrixMemoryPool_l), INTENT(INOUT) :: this
      INTEGER(kind=c_int), INTENT(IN) :: columns
      INTEGER(kind=c_int), INTENT(IN) :: rows
      REAL(NTREAL), INTENT(IN), OPTIONAL :: sparsity_in
    END SUBROUTINE Construct_base
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Destruct_base(this)
      USE DataTypesModule, ONLY : NTREAL
      IMPORT :: MatrixMemoryPool_l
      IMPLICIT NONE
      CLASS(MatrixMemoryPool_l), INTENT(INOUT) :: this
    END SUBROUTINE Destruct_base
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    PURE FUNCTION CheckValidity_base(this, columns, rows) RESULT(isvalid)
      IMPORT :: MatrixMemoryPool_l
      IMPLICIT NONE
      CLASS(MatrixMemoryPool_l), INTENT(IN) :: this
      INTEGER, INTENT(IN) :: columns
      INTEGER, INTENT(IN) :: rows
      LOGICAL :: isvalid
    END FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE SetSparsity_base(this, sparsity)
      USE DataTypesModule, ONLY : NTREAL
      IMPORT :: MatrixMemoryPool_l
      IMPLICIT NONE
      CLASS(MatrixMemoryPool_l), INTENT(INOUT) :: this
      REAL(NTREAL), INTENT(IN) :: sparsity
    END SUBROUTINE
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct Matrix Memory Pool object.
  !> @param[out] this a constructed Matrix Memory Pool object.
  !> @param[in] columns number of columns in the matrix.
  !> @param[in] rows number of rows in the matrix.
  !> @param[in] sparsity_in estimated sparsity (optional).
  SUBROUTINE Construct_r(this, columns, rows, sparsity_in)
    CLASS(MatrixMemoryPool_lr), INTENT(INOUT) :: this
    INTEGER(kind=c_int), INTENT(IN) :: columns
    INTEGER(kind=c_int), INTENT(IN) :: rows
    REAL(NTREAL), INTENT(IN), OPTIONAL :: sparsity_in

    INCLUDE "dense_includes/ConstructMatrixMemoryPool.f90"
  END SUBROUTINE Construct_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct Matrix Memory Pool object.
  !> @param[out] this a constructed Matrix Memory Pool object.
  !> @param[in] columns number of columns in the matrix.
  !> @param[in] rows number of rows in the matrix.
  !> @param[in] sparsity_in estimated sparsity (optional).
  SUBROUTINE Construct_c(this, columns, rows, sparsity_in)
    CLASS(MatrixMemoryPool_lc), INTENT(INOUT) :: this
    INTEGER(kind=c_int), INTENT(IN) :: columns
    INTEGER(kind=c_int), INTENT(IN) :: rows
    REAL(NTREAL), INTENT(IN), OPTIONAL :: sparsity_in

    INCLUDE "dense_includes/ConstructMatrixMemoryPool.f90"
  END SUBROUTINE Construct_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A destructor for a matrix memory pool
  !> @param[inout] this the matrix being destructed.
  SUBROUTINE Destruct_r(this)
    CLASS(MatrixMemoryPool_lr), INTENT(INOUT) :: this

    INCLUDE "dense_includes/DestructMatrixMemoryPool.f90"
  END SUBROUTINE Destruct_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A destructor for a matrix memory pool
  !> @param[inout] this the matrix being destructed.
  SUBROUTINE Destruct_c(this)
    CLASS(MatrixMemoryPool_lc), INTENT(INOUT) :: this

    INCLUDE "dense_includes/DestructMatrixMemoryPool.f90"
  END SUBROUTINE Destruct_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Checks if a given memory pool has been validly allocated to handle
  !! the given parameters.
  !> @param[in] this the memory pool to check.
  !> @param[in] columns number of columns in the matrix.
  !> @param[in] rows number of rows in the matrix.
  !> @return true if the memory pool is valid.
  PURE FUNCTION CheckValidity_r(this, columns, rows) RESULT(isvalid)
    CLASS(MatrixMemoryPool_lr), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: columns
    INTEGER, INTENT(IN) :: rows
    LOGICAL :: isvalid

    INCLUDE "dense_includes/CheckMemoryPoolValidity.f90"
  END FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Checks if a given memory pool has been validly allocated to handle
  !! the given parameters.
  !> @param[in] this the memory pool to check.
  !> @param[in] columns number of columns in the matrix.
  !> @param[in] rows number of rows in the matrix.
  !> @return true if the memory pool is valid.
  PURE FUNCTION CheckValidity_c(this, columns, rows) RESULT(isvalid)
    CLASS(MatrixMemoryPool_lc), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: columns
    INTEGER, INTENT(IN) :: rows
    LOGICAL :: isvalid

    INCLUDE "dense_includes/CheckMemoryPoolValidity.f90"
  END FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets the expected sparsity of the matrix, which helps with hashing.
  !! @param[inout] this the memory pool to set the sparsity of.
  !! @param[in] sparsity the sparsity value.
  SUBROUTINE SetSparsity_r(this, sparsity)
    CLASS(MatrixMemoryPool_lr), INTENT(INOUT) :: this
    REAL(NTREAL), INTENT(IN) :: sparsity

    INCLUDE "dense_includes/SetPoolSparsity.f90"
  END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets the expected sparsity of the matrix, which helps with hashing.
  !! @param[inout] this the memory pool to set the sparsity of.
  !! @param[in] sparsity the sparsity value.
  SUBROUTINE SetSparsity_c(this, sparsity)
    CLASS(MatrixMemoryPool_lc), INTENT(INOUT) :: this
    REAL(NTREAL), INTENT(IN) :: sparsity

    INCLUDE "dense_includes/SetPoolSparsity.f90"
  END SUBROUTINE
END MODULE MatrixMemoryPoolModule
