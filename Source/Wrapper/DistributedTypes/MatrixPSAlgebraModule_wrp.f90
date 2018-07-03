!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for wrapping a Distributed Sparse Matrix.
MODULE MatrixPSAlgebraModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE MatrixPSModule_wrp, ONLY : Matrix_ps_wrp
  USE MatrixMemoryPoolPModule_wrp, ONLY : MatrixMemoryPool_p_wrp
  USE MatrixPSAlgebraModule
  USE MatrixPSModule
  USE PermutationModule_wrp, ONLY : Permutation_wrp
  USE TripletListModule_wrp, ONLY : TripletList_r_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int, c_char, c_bool
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: IncrementMatrix_ps_wrp
  PUBLIC :: DotMatrix_ps_wrp
  PUBLIC :: MatrixPairwiseMultiply_ps_wrp
  PUBLIC :: MatrixMultiply_ps_wrp
  PUBLIC :: ScaleMatrix_ps_wrp
  PUBLIC :: MatrixNorm_ps_wrp
  PUBLIC :: MatrixTrace_ps_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Matrix B = alpha*Matrix A + Matrix B (AXPY)
  SUBROUTINE IncrementMatrix_ps_wrp(ih_matA, ih_matB, alpha_in,threshold_in) &
       & bind(c,name="IncrementMatrix_ps_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matB(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: alpha_in
    REAL(NTREAL), INTENT(IN) :: threshold_in
    TYPE(Matrix_ps_wrp) :: h_matA
    TYPE(Matrix_ps_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    CALL IncrementMatrix(h_matA%data, h_matB%data, alpha_in, threshold_in)
  END SUBROUTINE IncrementMatrix_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> product = dot(matA,matB)
  FUNCTION DotMatrix_ps_wrp(ih_matA, ih_matB) RESULT(product) &
       & bind(c,name="DotMatrix_ps_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matB(SIZE_wrp)
    REAL(NTREAL) :: product
    TYPE(Matrix_ps_wrp) :: h_matA
    TYPE(Matrix_ps_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    CALL DotMatrix(h_matA%data, h_matB%data, product)
  END FUNCTION DotMatrix_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Elementwise multiplication.
  SUBROUTINE MatrixPairwiseMultiply_ps_wrp(ih_matA, ih_matB, ih_matC) &
       & bind(c,name="MatrixPairwiseMultiply_ps_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_matB(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matC(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_matA
    TYPE(Matrix_ps_wrp) :: h_matB
    TYPE(Matrix_ps_wrp) :: h_matC

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    h_matC = TRANSFER(ih_matC,h_matC)

    CALL PairwiseMultiplyMatrix(h_matA%data, h_matB%data, h_matC%data)
  END SUBROUTINE MatrixPairwiseMultiply_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Multiply two matrices together, and add to the third.
  SUBROUTINE MatrixMultiply_ps_wrp(ih_matA, ih_matB, ih_matC, alpha_in, &
       & beta_in, threshold_in, ih_memory_pool_in) &
       & bind(c,name="MatrixMultiply_ps_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_matB(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matC(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: alpha_in
    REAL(NTREAL), INTENT(IN) :: beta_in
    REAL(NTREAL), INTENT(IN) :: threshold_in
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_memory_pool_in(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_matA
    TYPE(Matrix_ps_wrp) :: h_matB
    TYPE(Matrix_ps_wrp) :: h_matC
    TYPE(MatrixMemoryPool_p_wrp) :: h_memory_pool_in

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    h_matC = TRANSFER(ih_matC,h_matC)
    h_memory_pool_in = TRANSFER(ih_memory_pool_in,h_memory_pool_in)

    CALL MatrixMultiply(h_matA%data, h_matB%data, h_matC%data, &
         & alpha_in, beta_in, threshold_in, h_memory_pool_in%data)
  END SUBROUTINE MatrixMultiply_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a distributed sparse matrix by a constant.
  SUBROUTINE ScaleMatrix_ps_wrp(ih_this, constant) &
       & bind(c,name="ScaleMatrix_ps_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: constant
    TYPE(Matrix_ps_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL ScaleMatrix(h_this%data,constant)
  END SUBROUTINE ScaleMatrix_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the norm of a distributed sparse matrix along the rows.
  FUNCTION MatrixNorm_ps_wrp(ih_this) bind(c,name="MatrixNorm_ps_wrp") &
       & RESULT(norm_value)
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    REAL(NTREAL) :: norm_value
    TYPE(Matrix_ps_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    norm_value = MatrixNorm(h_this%data)
  END FUNCTION MatrixNorm_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the trace of the matrix.
  SUBROUTINE MatrixTrace_ps_wrp(ih_this, trace_value) &
        & bind(c,name="MatrixTrace_ps_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(OUT) :: trace_value
    TYPE(Matrix_ps_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL MatrixTrace(h_this%data, trace_value)
  END SUBROUTINE MatrixTrace_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE MatrixPSAlgebraModule_wrp
