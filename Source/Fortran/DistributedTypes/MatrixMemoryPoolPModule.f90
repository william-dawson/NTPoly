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
  PUBLIC :: ConstructMatrixMemoryPoolD
  PUBLIC :: DestructMatrixMemoryPoolD
  PUBLIC :: CheckMemoryPoolDValidity
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct Distributed Matrix Memory Pool object.
  !> @param[out] this a constructed Matrix Memory Pool object.
  !! @param[in] matrix the associated distributed sparse matrix.
  PURE SUBROUTINE ConstructMatrixMemoryPoolD(this, matrix)
    !! Parameters
    TYPE(MatrixMemoryPool_p), INTENT(INOUT) :: this
    TYPE(Matrix_ps), INTENT(IN) :: matrix

    !! Allocate
    CALL DestructMatrixMemoryPoolD(this)
    ALLOCATE(this%grid(matrix%process_grid%number_of_blocks_rows, &
         & matrix%process_grid%number_of_blocks_columns))
  END SUBROUTINE ConstructMatrixMemoryPoolD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a Distributed Matrix Memory Pool object.
  !> @param[out] this Distributed Matrix Memory Pool object to destroy.
  PURE SUBROUTINE DestructMatrixMemoryPoolD(this)
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
  END SUBROUTINE DestructMatrixMemoryPoolD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Checks if a given distributed memory pool has been validly allocated to
  !! handle the given parameters.
  !> @param[in] this the memory pool to check.
  !> @return true if the memory pool is valid.
  PURE FUNCTION CheckMemoryPoolDValidity(this) RESULT(isvalid)
    !! Parameters
    TYPE(MatrixMemoryPool_p), INTENT(IN) :: this
    LOGICAL :: isvalid

    isvalid = .TRUE.
    !! Check allocation
    IF (.NOT. ALLOCATED(this%grid)) isvalid = .FALSE.

  END FUNCTION CheckMemoryPoolDValidity
END MODULE MatrixMemoryPoolPModule
