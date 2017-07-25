!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for wrapping the distributed sparse matrix memory pool.
MODULE DistributedMatrixMemoryPoolModule_wrp
  USE DistributedMatrixMemoryPoolModule, ONLY : DistributedMatrixMemoryPool_t, &
       & ConstructDistributedMatrixMemoryPool, &
       & DestructDistributedMatrixMemoryPool
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrapper for the distributed matrix memory pool type
  TYPE, PUBLIC :: DistributedMatrixMemoryPool_wrp
     !> Actual data.
     TYPE(DistributedMatrixMemoryPool_t), POINTER :: data
  END TYPE DistributedMatrixMemoryPool_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructDistributedMatrixMemoryPool_wrp
  PUBLIC :: DestructDistributedMatrixMemoryPool_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct Distributed Matrix Memory Pool object.
  !> @param[out] ih_this to a constructed Matrix Memory Pool object.
  PURE SUBROUTINE ConstructDistributedMatrixMemoryPool_wrp(ih_this) &
       & bind(c,name="ConstructDistributedMatrixMemoryPool_wrp")
    !! Parameters
    INTEGER(kind=c_int), INTENT(out) :: ih_this(SIZE_wrp)
    TYPE(DistributedMatrixMemoryPool_wrp) :: h_this

    ALLOCATE(h_this%data)
    CALL ConstructDistributedMatrixMemoryPool(h_this%data)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructDistributedMatrixMemoryPool_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a permutation object.
  !! @param[inout] ih_this handle to the permutation to destruct.
  PURE SUBROUTINE DestructDistributedMatrixMemoryPool_wrp(ih_this) &
       & bind(c,name="DestructDistributedMatrixMemoryPool_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    TYPE(DistributedMatrixMemoryPool_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL DestructDistributedMatrixMemoryPool(h_this%data)
  END SUBROUTINE DestructDistributedMatrixMemoryPool_wrp
END MODULE DistributedMatrixMemoryPoolModule_wrp
