!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for wrapping the distributed sparse matrix memory pool.
MODULE PMatrixMemoryPoolModule_wrp
  USE PSMatrixModule_wrp, ONLY : Matrix_ps_wrp
  USE PMatrixMemoryPoolModule
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrapper for the distributed matrix memory pool type
  TYPE, PUBLIC :: MatrixMemoryPool_p_wrp
     TYPE(MatrixMemoryPool_p), POINTER :: DATA
  END TYPE MatrixMemoryPool_p_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructMatrixMemoryPool_p_wrp
  PUBLIC :: DestructMatrixMemoryPool_p_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct Distributed Matrix Memory Pool object.
  SUBROUTINE ConstructMatrixMemoryPool_p_wrp(ih_this, ih_matrix) &
       & BIND(c,name="ConstructMatrixMemoryPool_p_wrp")
    !! Parameters
    INTEGER(kind=c_int), INTENT(out) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: ih_matrix(SIZE_wrp)
    TYPE(MatrixMemoryPool_p_wrp) :: h_this
    TYPE(Matrix_ps_wrp) :: h_matrix

    h_matrix = TRANSFER(ih_matrix,h_matrix)

    ALLOCATE(h_this%DATA)
    CALL ConstructMatrixMemoryPool(h_this%DATA, h_matrix%DATA)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructMatrixMemoryPool_p_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a permutation object.
  SUBROUTINE DestructMatrixMemoryPool_p_wrp(ih_this) &
       & BIND(c,name="DestructMatrixMemoryPool_p_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    TYPE(MatrixMemoryPool_p_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL DestructMatrixMemoryPool(h_this%DATA)
  END SUBROUTINE DestructMatrixMemoryPool_p_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE PMatrixMemoryPoolModule_wrp
