!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling scratch memory for distributed matrix multiplication.
MODULE DistributedMatrixMemoryPoolModule
  USE ErrorModule
  USE MatrixMemoryPoolModule, ONLY : MatrixMemoryPool_t, &
       & DestructMatrixMemoryPool
  USE ProcessGridModule, ONLY : number_of_blocks_columns, number_of_blocks_rows
  USE iso_c_binding, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A memory pool datatype that can be reused for matrix matrix multiplication.
  !! this is to prevent excessive alloc/dealloc.
  TYPE, PUBLIC :: DistributedMatrixMemoryPool_t
     !> Grid of local pools.
     TYPE(MatrixMemoryPool_t), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: grid
  END TYPE DistributedMatrixMemoryPool_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructDistributedMatrixMemoryPool
  PUBLIC :: DestructDistributedMatrixMemoryPool
  PUBLIC :: CheckDistributedMemoryPoolValidity
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct Distributed Matrix Memory Pool object.
  !> @param[out] this a constructed Matrix Memory Pool object.
  PURE SUBROUTINE ConstructDistributedMatrixMemoryPool(this)
    !! Parameters
    TYPE(DistributedMatrixMemoryPool_t), INTENT(out) :: this

    !! Allocate
    CALL DestructDistributedMatrixMemoryPool(this)
    ALLOCATE(this%grid(number_of_blocks_columns,number_of_blocks_rows))
  END SUBROUTINE ConstructDistributedMatrixMemoryPool
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a Distributed Matrix Memory Pool object.
  !> @param[out] this Distributed Matrix Memory Pool object to destroy.
  PURE SUBROUTINE DestructDistributedMatrixMemoryPool(this)
    !! Parameters
    TYPE(DistributedMatrixMemoryPool_t), INTENT(out) :: this
    !! Local Data
    INTEGER :: row_counter, column_counter

    !! Allocate
    IF (ALLOCATED(this%grid)) THEN
       DO row_counter = LBOUND(this%grid,2), UBOUND(this%grid,2)
          DO column_counter = LBOUND(this%grid,1), UBOUND(this%grid,1)
             CALL DestructMatrixMemoryPool( &
                  & this%grid(column_counter,row_counter))
          END DO
       END DO
       DEALLOCATE(this%grid)
    END IF
  END SUBROUTINE DestructDistributedMatrixMemoryPool
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Checks if a given distributed memory pool has been validly allocated to handle
  !! the given parameters.
  !> @param[in] this the memory pool to check.
  !> @return true if the memory pool is valid.
  PURE FUNCTION CheckDistributedMemoryPoolValidity(this) &
       & RESULT(isvalid)
    !! Parameters
    TYPE(DistributedMatrixMemoryPool_t), INTENT(in) :: this
    LOGICAL :: isvalid

    isvalid = .TRUE.
    !! Check allocation
    IF (.NOT. ALLOCATED(this%grid)) isvalid = .FALSE.

  END FUNCTION CheckDistributedMemoryPoolValidity
END MODULE DistributedMatrixMemoryPoolModule
