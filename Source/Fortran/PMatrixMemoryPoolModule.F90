!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling scratch memory for distributed matrix multiplication.
MODULE PMatrixMemoryPoolModule
#ifdef NTPOLY_ENABLE_GPU_TRS2
  USE ISO_C_BINDING, ONLY: C_PTR,C_NULL_PTR,C_ASSOCIATED,C_INT
#endif
  USE PSMatrixModule, ONLY : Matrix_ps
  USE MatrixMemoryPoolModule, ONLY : MatrixMemoryPool_lr, MatrixMemoryPool_lc, &
       & DestructMatrixMemoryPool
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A memory pool datatype that can be reused for matrix matrix multiplication.
  !> this is to prevent excessive alloc/dealloc.
  TYPE, PUBLIC :: MatrixMemoryPool_p
     !> Grid of local pools.
     TYPE(MatrixMemoryPool_lr), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: grid_r
     !> Grid of local pools (complex).
     TYPE(MatrixMemoryPool_lc), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: grid_c
#ifdef NTPOLY_ENABLE_GPU_TRS2
     TYPE(C_PTR), PRIVATE :: gpu_workspace=C_NULL_PTR
   CONTAINS
     FINAL :: FinalizeMemoryPool
     PROCEDURE, PRIVATE :: AssignMemoryPool
     GENERIC :: ASSIGNMENT(=) => AssignMemoryPool
#endif
  END TYPE MatrixMemoryPool_p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructMatrixMemoryPool
  PUBLIC :: DestructMatrixMemoryPool
  PUBLIC :: CheckMemoryPoolValidity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE ConstructMatrixMemoryPool
     MODULE PROCEDURE ConstructMatrixMemoryPool_p
  END INTERFACE ConstructMatrixMemoryPool
  INTERFACE DestructMatrixMemoryPool
     MODULE PROCEDURE DestructMatrixMemoryPool_p
  END INTERFACE DestructMatrixMemoryPool
  INTERFACE CheckMemoryPoolValidity
     MODULE PROCEDURE CheckMemoryPoolValidity_p
  END INTERFACE CheckMemoryPoolValidity
#ifdef NTPOLY_ENABLE_GPU_TRS2
  PUBLIC :: ActivateGPUWorkspace
  INTERFACE
    INTEGER(C_INT) FUNCTION select_gpu_pool(owner) BIND(C,NAME='solver_gpu_pool_select_')
      IMPORT C_PTR,C_INT
      TYPE(C_PTR) :: owner
    END FUNCTION
    INTEGER(C_INT) FUNCTION retain_gpu_pool(owner) BIND(C,NAME='solver_gpu_pool_retain_')
      IMPORT C_PTR,C_INT
      TYPE(C_PTR), VALUE :: owner
    END FUNCTION
    INTEGER(C_INT) FUNCTION release_gpu_pool(owner) BIND(C,NAME='solver_gpu_pool_release_')
      IMPORT C_PTR,C_INT
      TYPE(C_PTR) :: owner
    END FUNCTION
  END INTERFACE
#endif
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct Distributed Matrix Memory Pool object.
#ifdef NTPOLY_ENABLE_GPU_TRS2
  SUBROUTINE ConstructMatrixMemoryPool_p(this, matrix)
#else
  PURE SUBROUTINE ConstructMatrixMemoryPool_p(this, matrix)
#endif
    !> A constructed Matrix Memory Pool object.
    TYPE(MatrixMemoryPool_p), INTENT(INOUT) :: this
    !> The associated distributed sparse matrix.
    TYPE(Matrix_ps), INTENT(IN) :: matrix

    CALL DestructMatrixMemoryPool(this)
    !! Allocate
    IF (matrix%is_complex) THEN
       ALLOCATE(this%grid_c(matrix%process_grid%number_of_blocks_rows, &
            & matrix%process_grid%number_of_blocks_columns))
    ELSE
       ALLOCATE(this%grid_r(matrix%process_grid%number_of_blocks_rows, &
            & matrix%process_grid%number_of_blocks_columns))
    END IF
  END SUBROUTINE ConstructMatrixMemoryPool_p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a Distributed Matrix Memory Pool object.
#ifdef NTPOLY_ENABLE_GPU_TRS2
  SUBROUTINE DestructMatrixMemoryPool_p(this)
#else
  PURE SUBROUTINE DestructMatrixMemoryPool_p(this)
#endif
    !> Distributed Matrix Memory Pool object to destroy.
    CLASS(MatrixMemoryPool_p), INTENT(INOUT) :: this
    !! Local Data
    INTEGER :: II, JJ

#define grid grid_r
#include "distributed_pool_includes/DestructMatrixMemoryPool.f90"
#undef grid

#define grid grid_c
#include "distributed_pool_includes/DestructMatrixMemoryPool.f90"
#undef grid

#ifdef NTPOLY_ENABLE_GPU_TRS2
    IF(C_ASSOCIATED(this%gpu_workspace)) THEN
       IF(release_gpu_pool(this%gpu_workspace)/=0) ERROR STOP 'GPU memory pool release failed'
    END IF
#endif
  END SUBROUTINE DestructMatrixMemoryPool_p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Checks if a given distributed memory pool has been validly allocated to
  !> handle the given parameters.
  PURE FUNCTION CheckMemoryPoolValidity_p(this, matrix) RESULT(isvalid)
    !> The memory pool to check.
    TYPE(MatrixMemoryPool_p), INTENT(IN) :: this
    !> The associated matrix to check against.
    TYPE(Matrix_ps), INTENT(IN) :: matrix
    !> True if the memory pool is valid.
    LOGICAL :: isvalid

    isvalid = .TRUE.
    !! Check allocation
    IF (matrix%is_complex) THEN
       IF (.NOT. ALLOCATED(this%grid_c)) isvalid = .FALSE.
    ELSE
       IF (.NOT. ALLOCATED(this%grid_r)) isvalid = .FALSE.
    END IF

  END FUNCTION CheckMemoryPoolValidity_p
#ifdef NTPOLY_ENABLE_GPU_TRS2
  ! Scratch storage is shared on assignment, with independent CPU pool arrays.
  ! GPU calls are serialized by the context lease and refresh input payloads.
  IMPURE ELEMENTAL SUBROUTINE AssignMemoryPool(this,other)
    CLASS(MatrixMemoryPool_p), INTENT(INOUT) :: this
    TYPE(MatrixMemoryPool_p), INTENT(IN) :: other
    TYPE(MatrixMemoryPool_lr), ALLOCATABLE :: real_grid(:,:)
    TYPE(MatrixMemoryPool_lc), ALLOCATABLE :: complex_grid(:,:)
    TYPE(C_PTR) :: owner
    IF(ALLOCATED(other%grid_r)) real_grid=other%grid_r
    IF(ALLOCATED(other%grid_c)) complex_grid=other%grid_c
    owner=other%gpu_workspace
    IF(C_ASSOCIATED(owner)) THEN
       IF(retain_gpu_pool(owner)/=0) ERROR STOP 'GPU memory pool retain failed'
    END IF
    CALL DestructMatrixMemoryPool_p(this)
    CALL MOVE_ALLOC(real_grid,this%grid_r)
    CALL MOVE_ALLOC(complex_grid,this%grid_c)
    this%gpu_workspace=owner
  END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE FinalizeMemoryPool(this)
    TYPE(MatrixMemoryPool_p), INTENT(INOUT) :: this
    CALL DestructMatrixMemoryPool_p(this)
  END SUBROUTINE

  SUBROUTINE ActivateGPUWorkspace(this)
    TYPE(MatrixMemoryPool_p), INTENT(INOUT) :: this
    IF(select_gpu_pool(this%gpu_workspace)/=0) ERROR STOP 'GPU memory pool activation failed'
  END SUBROUTINE
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE PMatrixMemoryPoolModule
