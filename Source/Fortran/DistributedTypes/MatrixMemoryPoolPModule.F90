!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling scratch memory for distributed matrix multiplication.
MODULE MatrixMemoryPoolPModule
  USE MatrixPSModule, ONLY : Matrix_ps
  USE MatrixMemoryPoolModule, ONLY : MatrixMemoryPool_lr, &
       & DestructMatrixMemoryPool
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A memory pool datatype that can be reused for matrix matrix multiplication.
  !! this is to prevent excessive alloc/dealloc.
  TYPE, PUBLIC :: MatrixMemoryPool_p
     !> Grid of local pools.
     TYPE(MatrixMemoryPool_lr), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: grid
  END TYPE MatrixMemoryPool_p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: DestructMatrixMemoryPool
  PUBLIC :: CheckMemoryPoolValidity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE MatrixMemoryPool_p
    MODULE PROCEDURE ConstructMatrixMemoryPool_p
  END INTERFACE
  INTERFACE DestructMatrixMemoryPool
    MODULE PROCEDURE DestructMatrixMemoryPool_p
  END INTERFACE
  INTERFACE CheckMemoryPoolValidity
    MODULE PROCEDURE CheckMemoryPoolValidity_p
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct Distributed Matrix Memory Pool object.
  !> @param[out] this a constructed Matrix Memory Pool object.
  !! @param[in] matrix the associated distributed sparse matrix.
  PURE FUNCTION ConstructMatrixMemoryPool_p(matrix) RESULT(this)
    !! Parameters
    TYPE(MatrixMemoryPool_p) :: this
    TYPE(Matrix_ps), INTENT(IN) :: matrix

    !! Allocate
    ALLOCATE(this%grid(matrix%process_grid%number_of_blocks_rows, &
         & matrix%process_grid%number_of_blocks_columns))
  END FUNCTION ConstructMatrixMemoryPool_p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a Distributed Matrix Memory Pool object.
  !> @param[out] this Distributed Matrix Memory Pool object to destroy.
  PURE SUBROUTINE DestructMatrixMemoryPool_p(this)
    !! Parameters
    TYPE(MatrixMemoryPool_p), INTENT(INOUT) :: this
    !! Local Data
    INTEGER :: row_counter, column_counter

    !! Allocate
    IF (ALLOCATED(this%grid)) THEN
       DO column_counter = LBOUND(this%grid,2), UBOUND(this%grid,2)
          DO row_counter = LBOUND(this%grid,1), UBOUND(this%grid,1)
             CALL DestructMatrixMemoryPool( &
                  & this%grid(row_counter, column_counter))
          END DO
       END DO
       DEALLOCATE(this%grid)
    END IF
  END SUBROUTINE DestructMatrixMemoryPool_p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Checks if a given distributed memory pool has been validly allocated to
  !! handle the given parameters.
  !> @param[in] this the memory pool to check.
  !> @return true if the memory pool is valid.
  PURE FUNCTION CheckMemoryPoolValidity_p(this) RESULT(isvalid)
    !! Parameters
    TYPE(MatrixMemoryPool_p), INTENT(IN) :: this
    LOGICAL :: isvalid

    isvalid = .TRUE.
    !! Check allocation
    IF (.NOT. ALLOCATED(this%grid)) isvalid = .FALSE.

  END FUNCTION CheckMemoryPoolValidity_p
END MODULE MatrixMemoryPoolPModule
