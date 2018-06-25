!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling scratch memory for distributed matrix multiplication.
MODULE MatrixMemoryPoolPModule
  USE MatrixPSModule, ONLY : Matrix_ps, Matrix_psr
  USE MatrixMemoryPoolModule, ONLY : MatrixMemoryPool_lr, MatrixMemoryPool_lc
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A memory pool datatype that can be reused for matrix matrix multiplication.
  !! this is to prevent excessive alloc/dealloc.
  TYPE, ABSTRACT, PUBLIC :: MatrixMemoryPool_p
   CONTAINS
     PROCEDURE(Construct_base), DEFERRED :: Init
     PROCEDURE(Destruct_base), DEFERRED :: Destruct
     PROCEDURE(CheckValidity_base), DEFERRED :: CheckValidity
  END TYPE MatrixMemoryPool_p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE, EXTENDS(MatrixMemoryPool_p), PUBLIC :: MatrixMemoryPool_pr
     !> Grid of local pools.
     TYPE(MatrixMemoryPool_lr), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: grid
   CONTAINS
     PROCEDURE :: Init => Construct_pr
     PROCEDURE :: Destruct => Destruct_pr
     PROCEDURE :: CheckValidity => CheckValidity_pr
  END TYPE MatrixMemoryPool_pr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! TYPE, EXTENDS(MatrixMemoryPool_p), PUBLIC :: MatrixMemoryPool_pc
  !    !> Grid of local pools.
  !    TYPE(MatrixMemoryPool_lc), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: grid
  !  CONTAINS
  !    PROCEDURE :: Init => Construct_pc
  !    PROCEDURE :: Destruct => Destruct_pc
  !    PROCEDURE :: CheckValidity => CheckValidity_pc
  ! END TYPE MatrixMemoryPool_pc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ABSTRACT INTERFACE
    PURE SUBROUTINE Construct_base(this, matrix)
      USE MatrixPSModule, ONLY : Matrix_ps
      IMPORT :: MatrixMemoryPool_p
      IMPLICIT NONE
      CLASS(MatrixMemoryPool_p), INTENT(INOUT) :: this
      CLASS(Matrix_ps), INTENT(IN) :: matrix
    END SUBROUTINE Construct_base
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    PURE SUBROUTINE Destruct_base(this)
      USE MatrixPSModule, ONLY : Matrix_ps
      IMPORT :: MatrixMemoryPool_p
      IMPLICIT NONE
      CLASS(MatrixMemoryPool_p), INTENT(INOUT) :: this
    END SUBROUTINE Destruct_base
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    PURE FUNCTION CheckValidity_base(this) RESULT(isvalid)
      USE MatrixPSModule, ONLY : Matrix_ps
      IMPORT :: MatrixMemoryPool_p
      IMPLICIT NONE
      CLASS(MatrixMemoryPool_p), INTENT(IN) :: this
      LOGICAL :: isvalid
    END FUNCTION CheckValidity_base
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct Distributed Matrix Memory Pool object.
  !> @param[out] this a constructed Matrix Memory Pool object.
  !! @param[in] matrix the associated distributed sparse matrix.
  PURE SUBROUTINE Construct_pr(this, matrix)
    !! Parameters
    CLASS(MatrixMemoryPool_pr), INTENT(INOUT) :: this
    CLASS(Matrix_ps), INTENT(IN) :: matrix

    INCLUDE "includes/ConstructMatrixMemoryPool.f90"
  END SUBROUTINE Construct_pr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !> Construct Distributed Matrix Memory Pool object.
  ! !> @param[out] this a constructed Matrix Memory Pool object.
  ! !! @param[in] matrix the associated distributed sparse matrix.
  ! PURE SUBROUTINE Construct_pc(this, matrix)
  !   !! Parameters
  !   CLASS(MatrixMemoryPool_pc), INTENT(INOUT) :: this
  !   CLASS(Matrix_psc), INTENT(IN) :: matrix
  !
  !   INCLUDE "includes/ConstructMatrixMemoryPool.f90"
  ! END SUBROUTINE Construct_pc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a Distributed Matrix Memory Pool object.
  !> @param[out] this Distributed Matrix Memory Pool object to destroy.
  PURE SUBROUTINE Destruct_pr(this)
    !! Parameters
    CLASS(MatrixMemoryPool_pr), INTENT(INOUT) :: this
    !! Local Data
    INTEGER :: row_counter, column_counter

    INCLUDE "includes/DestructMatrixMemoryPool.f90"
  END SUBROUTINE Destruct_pr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !> Destruct a Distributed Matrix Memory Pool object.
  ! !> @param[out] this Distributed Matrix Memory Pool object to destroy.
  ! PURE SUBROUTINE Destruct_pc(this)
  !   !! Parameters
  !   CLASS(MatrixMemoryPool_pc), INTENT(INOUT) :: this
  !   !! Local Data
  !   INTEGER :: row_counter, column_counter
  !
  !   INCLUDE "includes/DestructMatrixMemoryPool.f90"
  ! END SUBROUTINE Destruct_pc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Checks if a given distributed memory pool has been validly allocated to
  !! handle the given parameters.
  !> @param[in] this the memory pool to check.
  !> @return true if the memory pool is valid.
  PURE FUNCTION CheckValidity_pr(this) RESULT(isvalid)
    !! Parameters
    CLASS(MatrixMemoryPool_pr), INTENT(IN) :: this
    LOGICAL :: isvalid

    INCLUDE "includes/CheckValidity.f90"
  END FUNCTION CheckValidity_pr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !> Checks if a given distributed memory pool has been validly allocated to
  ! !! handle the given parameters.
  ! !> @param[in] this the memory pool to check.
  ! !> @return true if the memory pool is valid.
  ! PURE FUNCTION CheckValidity_pc(this) RESULT(isvalid)
  !   !! Parameters
  !   TYPE(MatrixMemoryPool_pc), INTENT(IN) :: this
  !   LOGICAL :: isvalid
  !
  !   INCLUDE "includes/CheckValidity.f90"
  ! END FUNCTION CheckValidity_pc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE MatrixMemoryPoolPModule
