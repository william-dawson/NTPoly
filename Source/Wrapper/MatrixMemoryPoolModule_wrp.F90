!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for wrapping the matrix memory pool data type.
MODULE MatrixMemoryPoolModule_wrp
  USE MatrixMemoryPoolModule, ONLY : MatrixMemoryPool_lr, &
       & MatrixMemoryPool_lc, DestructMatrixMemoryPool, &
       & ConstructMatrixMemoryPool
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
  !> A wrapper for the matrix memory pool data type.
  TYPE, PUBLIC :: MatrixMemoryPool_lc_wrp
     TYPE(MatrixMemoryPool_lc), POINTER :: DATA
  END TYPE MatrixMemoryPool_lc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructMatrixMemoryPool_lr_wrp
  PUBLIC :: DestructMatrixMemoryPool_lr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructMatrixMemoryPool_lc_wrp
  PUBLIC :: DestructMatrixMemoryPool_lc_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the Matrix Memory Pool constructor.
  SUBROUTINE ConstructMatrixMemoryPool_lr_wrp(ih_this, columns, rows) &
       & BIND(c,name="ConstructMatrixMemoryPool_lr_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: columns
    INTEGER(kind=c_int), INTENT(in) :: rows
    TYPE(MatrixMemoryPool_lr_wrp) :: h_this

    ALLOCATE(h_this%DATA)
    CALL ConstructMatrixMemoryPool(h_this%DATA, columns, rows)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructMatrixMemoryPool_lr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the destructor for a matrix memory pool
  SUBROUTINE DestructMatrixMemoryPool_lr_wrp(ih_this) &
       & BIND(c,name="DestructMatrixMemoryPool_lr_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    TYPE(MatrixMemoryPool_lr_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL DestructMatrixMemoryPool(h_this%DATA)
    DEALLOCATE(h_this%DATA)
    !ih_this = 0
  END SUBROUTINE DestructMatrixMemoryPool_lr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the Matrix Memory Pool constructor.
  SUBROUTINE ConstructMatrixMemoryPool_lc_wrp(ih_this, columns, rows) &
       & BIND(c,name="ConstructMatrixMemoryPool_lc_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: columns
    INTEGER(kind=c_int), INTENT(in) :: rows
    TYPE(MatrixMemoryPool_lc_wrp) :: h_this

    ALLOCATE(h_this%DATA)
    CALL ConstructMatrixMemoryPool(h_this%DATA, columns, rows)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructMatrixMemoryPool_lc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the destructor for a matrix memory pool
  SUBROUTINE DestructMatrixMemoryPool_lc_wrp(ih_this) &
       & BIND(c,name="DestructMatrixMemoryPool_lc_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    TYPE(MatrixMemoryPool_lc_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL DestructMatrixMemoryPool(h_this%DATA)
    DEALLOCATE(h_this%DATA)
    !ih_this = 0
  END SUBROUTINE DestructMatrixMemoryPool_lc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE MatrixMemoryPoolModule_wrp
