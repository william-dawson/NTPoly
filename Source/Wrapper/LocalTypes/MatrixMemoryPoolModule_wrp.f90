!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for wrapping the matrix memory pool data type.
MODULE MatrixMemoryPoolModule_wrp
  USE MatrixMemoryPoolModule, ONLY : MatrixMemoryPool_lr, &
       & ConstructMatrixMemoryPool, DestructMatrixMemoryPool
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for the matrix memory pool data type.
  TYPE, PUBLIC :: MatrixMemoryPool_lr_wrp
     TYPE(MatrixMemoryPool_lr), POINTER :: DATA
  END TYPE MatrixMemoryPool_lr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructMatrixMemoryPool_lr_wrp
  PUBLIC :: DestructMatrixMemoryPool_lr_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the Matrix Memory Pool constructor.
  SUBROUTINE ConstructMatrixMemoryPool_lr_wrp(ih_this, columns, rows) &
       & bind(c,name="ConstructMatrixMemoryPool_lr_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: columns
    INTEGER(kind=c_int), INTENT(in) :: rows
    TYPE(MatrixMemoryPool_lr_wrp) :: h_this

    ALLOCATE(h_this%data)
    CALL ConstructMatrixMemoryPool(h_this%data,columns,rows)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructMatrixMemoryPool_lr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the destructor for a matrix memory pool
  PURE SUBROUTINE DestructMatrixMemoryPool_lr_wrp(ih_this) &
       & bind(c,name="DestructMatrixMemoryPool_lr_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    TYPE(MatrixMemoryPool_lr_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL DestructMatrixMemoryPool(h_this%data)
    DEALLOCATE(h_this%data)
    !ih_this = 0
  END SUBROUTINE DestructMatrixMemoryPool_lr_wrp
END MODULE MatrixMemoryPoolModule_wrp
