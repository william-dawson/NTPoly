!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for wrapping the sparse algebra routines.
MODULE MatrixSAlgebraModule_wrp
  USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
  USE MatrixMemoryPoolModule_wrp, ONLY : MatrixMemoryPool_lr_wrp, &
       & MatrixMemoryPool_lc_wrp
  USE MatrixSAlgebraModule, ONLY : ScaleMatrix, IncrementMatrix, &
       & MatrixMultiply, DotMatrix, PairwiseMultiplyMatrix
  USE MatrixSModule, ONLY : Matrix_lsr, Matrix_lsc
  USE MatrixSModule_wrp, ONLY: Matrix_lsr_wrp, Matrix_lsc_wrp
  USE TripletListModule_wrp, ONLY : TripletList_r_wrp, TripletList_c_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int, c_char, c_bool
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ScaleMatrix_lsr_wrp
  PUBLIC :: IncrementMatrix_lsr_wrp
  PUBLIC :: DotMatrix_lsr_wrp
  PUBLIC :: PairwiseMultiplyMatrix_lsr_wrp
  PUBLIC :: MatrixMultiply_lsr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ScaleMatrix_lsc_wrp
  PUBLIC :: IncrementMatrix_lsc_wrp
  PUBLIC :: DotMatrix_lsc_wrp
  PUBLIC :: PairwiseMultiplyMatrix_lsc_wrp
  PUBLIC :: MatrixMultiply_lsc_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the scale a sparse matrix by a constant routine.
  PURE SUBROUTINE ScaleMatrix_lsr_wrp(ih_this, constant) &
       & bind(c,name="ScaleMatrix_lsr_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: constant
    TYPE(Matrix_lsr_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL ScaleMatrix(h_this%data,constant)
  END SUBROUTINE ScaleMatrix_lsr_wrp
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap matrix incrementing function.
  PURE SUBROUTINE IncrementMatrix_lsr_wrp(ih_matA, ih_matB, alpha_in, &
       & threshold_in) bind(c,name="IncrementMatrix_lsr_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matB(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: alpha_in
    REAL(NTREAL), INTENT(IN) :: threshold_in
    TYPE(Matrix_lsr_wrp) :: h_matA
    TYPE(Matrix_lsr_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    CALL IncrementMatrix(h_matA%data, h_matB%data, alpha_in, threshold_in)
  END SUBROUTINE IncrementMatrix_lsr_wrp
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap matrix dot product function.
  PURE FUNCTION DotMatrix_lsr_wrp(ih_matA, ih_matB) RESULT(product) &
       & bind(c,name="DotMatrix_lsr_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_matB(SIZE_wrp)
    REAL(NTREAL) :: product
    TYPE(Matrix_lsr_wrp) :: h_matA
    TYPE(Matrix_lsr_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    product = DotMatrix(h_matA%data, h_matB%data)
  END FUNCTION DotMatrix_lsr_wrp
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap pairwise matrix multiplication function.
  SUBROUTINE PairwiseMultiplyMatrix_lsr_wrp(ih_matA, ih_matB, ih_matC) &
       & bind(c,name="PairwiseMultiplyMatrix_lsr_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_matB(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matC(SIZE_wrp)
    TYPE(Matrix_lsr_wrp) :: h_matA
    TYPE(Matrix_lsr_wrp) :: h_matB
    TYPE(Matrix_lsr_wrp) :: h_matC

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    h_matC = TRANSFER(ih_matC,h_matC)

    CALL PairwiseMultiplyMatrix(h_matA%data, h_matB%data, h_matC%data)
  END SUBROUTINE PairwiseMultiplyMatrix_lsr_wrp
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap matrix multiplication function.
  SUBROUTINE MatrixMultiply_lsr_wrp(ih_matA, ih_matB, ih_matC, IsATransposed, &
       & IsBTransposed, alpha, beta, threshold, ih_blocked_memory_pool) &
       & bind(c,name="MatrixMultiply_lsr_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_matB(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matC(SIZE_wrp)
    LOGICAL(kind=c_bool), INTENT(IN) :: IsATransposed
    LOGICAL(kind=c_bool), INTENT(IN) :: IsBTransposed
    REAL(NTREAL), INTENT(in) :: alpha
    REAL(NTREAL), INTENT(in) :: beta
    REAL(NTREAL), INTENT(in) :: threshold
    INTEGER(kind=c_int), INTENT(inout) :: ih_blocked_memory_pool(SIZE_wrp)
    TYPE(Matrix_lsr_wrp) :: h_matA
    TYPE(Matrix_lsr_wrp) :: h_matB
    TYPE(Matrix_lsr_wrp) :: h_matC
    TYPE(MatrixMemoryPool_lr_wrp) :: h_blocked_memory_pool

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    h_matC = TRANSFER(ih_matC,h_matC)
    h_blocked_memory_pool = TRANSFER(ih_blocked_memory_pool, &
         & h_blocked_memory_pool)

    CALL MatrixMultiply(h_matA%data, h_matB%data, h_matC%data, &
         & LOGICAL(IsATransposed), LOGICAL(IsBTransposed), alpha, &
         & beta, threshold, h_blocked_memory_pool%data)
  END SUBROUTINE MatrixMultiply_lsr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the scale a sparse matrix by a constant routine.
  PURE SUBROUTINE ScaleMatrix_lsc_wrp(ih_this, constant) &
       & bind(c,name="ScaleMatrix_lsc_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: constant
    TYPE(Matrix_lsc_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL ScaleMatrix(h_this%data,constant)
  END SUBROUTINE ScaleMatrix_lsc_wrp
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap matrix incrementing function.
  PURE SUBROUTINE IncrementMatrix_lsc_wrp(ih_matA, ih_matB, alpha_in, &
       & threshold_in) bind(c,name="IncrementMatrix_lsc_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matB(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: alpha_in
    REAL(NTREAL), INTENT(IN) :: threshold_in
    TYPE(Matrix_lsc_wrp) :: h_matA
    TYPE(Matrix_lsc_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    CALL IncrementMatrix(h_matA%data, h_matB%data, alpha_in, threshold_in)
  END SUBROUTINE IncrementMatrix_lsc_wrp
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap matrix dot product function.
  PURE FUNCTION DotMatrix_lsc_wrp(ih_matA, ih_matB) RESULT(product) &
       & bind(c,name="DotMatrix_lsc_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_matB(SIZE_wrp)
    REAL(NTREAL) :: product
    TYPE(Matrix_lsc_wrp) :: h_matA
    TYPE(Matrix_lsc_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    product = DotMatrix(h_matA%data, h_matB%data)
  END FUNCTION DotMatrix_lsc_wrp
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap pairwise matrix multiplication function.
  SUBROUTINE PairwiseMultiplyMatrix_lsc_wrp(ih_matA, ih_matB, ih_matC) &
       & bind(c,name="PairwiseMultiplyMatrix_lsc_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_matB(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matC(SIZE_wrp)
    TYPE(Matrix_lsc_wrp) :: h_matA
    TYPE(Matrix_lsc_wrp) :: h_matB
    TYPE(Matrix_lsc_wrp) :: h_matC

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    h_matC = TRANSFER(ih_matC,h_matC)

    CALL PairwiseMultiplyMatrix(h_matA%data, h_matB%data, h_matC%data)
  END SUBROUTINE PairwiseMultiplyMatrix_lsc_wrp
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap matrix multiplication function.
  SUBROUTINE MatrixMultiply_lsc_wrp(ih_matA, ih_matB, ih_matC, IsATransposed, &
       & IsBTransposed, alpha, beta, threshold, ih_blocked_memory_pool) &
       & bind(c,name="MatrixMultiply_lsc_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_matB(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matC(SIZE_wrp)
    LOGICAL(kind=c_bool), INTENT(IN) :: IsATransposed
    LOGICAL(kind=c_bool), INTENT(IN) :: IsBTransposed
    REAL(NTREAL), INTENT(in) :: alpha
    REAL(NTREAL), INTENT(in) :: beta
    REAL(NTREAL), INTENT(in) :: threshold
    INTEGER(kind=c_int), INTENT(inout) :: ih_blocked_memory_pool(SIZE_wrp)
    TYPE(Matrix_lsc_wrp) :: h_matA
    TYPE(Matrix_lsc_wrp) :: h_matB
    TYPE(Matrix_lsc_wrp) :: h_matC
    TYPE(MatrixMemoryPool_lc_wrp) :: h_blocked_memory_pool

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    h_matC = TRANSFER(ih_matC,h_matC)
    h_blocked_memory_pool = TRANSFER(ih_blocked_memory_pool, &
         & h_blocked_memory_pool)

    CALL MatrixMultiply(h_matA%data, h_matB%data, h_matC%data, &
         & LOGICAL(IsATransposed), LOGICAL(IsBTransposed), alpha, &
         & beta, threshold, h_blocked_memory_pool%data)
  END SUBROUTINE MatrixMultiply_lsc_wrp
END MODULE MatrixSAlgebraModule_wrp
