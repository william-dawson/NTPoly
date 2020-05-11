!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This module contains helper routines for converting an NTPoly matrix
!> to data structures used in other programs.
MODULE MatrixConversionModule_wrp
  USE MatrixConversionModule
  USE PSMatrixModule_wrp, ONLY : Matrix_ps_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : C_INT
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SnapMatrixToSparsityPattern_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Some codes use a fixed sparsity pattern for a matrix instead of filtering
  !> small values. Using this routine, the matrix is filled to have the same
  !> pattern as the second matrix argument. Zeros of the sparsity pattern are
  !> left in, whereas values outside the sparsity are removed. This can 
  !> faciliate conversion between formats.
  SUBROUTINE SnapMatrixToSparsityPattern_wrp(ih_matA, ih_matB) &
       & BIND(C,NAME="SnapMatrixToSparsityPattern_wrp")
    INTEGER(KIND=C_INT), INTENT(INOUT) :: ih_matA(SIZE_wrp)
    INTEGER(KIND=C_INT), INTENT(IN) :: ih_matB(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_matA
    TYPE(Matrix_ps_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    CALL SnapMatrixToSparsityPattern(h_matA%data, h_matB%data)
  END SUBROUTINE SnapMatrixToSparsityPattern_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE