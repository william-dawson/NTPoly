!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for wrapping a Distributed Sparse Matrix.
MODULE DistributedSparseMatrixAlgebraModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE DistributedSparseMatrixModule_wrp, ONLY : DistributedSparseMatrix_wrp
  USE DistributedMatrixMemoryPoolModule_wrp, ONLY : &
       & DistributedMatrixMemoryPool_wrp
  USE DistributedSparseMatrixAlgebraModule
  USE DistributedSparseMatrixModule
  USE PermutationModule_wrp, ONLY : Permutation_wrp
  USE TripletListModule_wrp, ONLY : TripletList_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int, c_char, c_bool
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: IncrementDistributedSparseMatrix_wrp
  PUBLIC :: DotDistributedSparseMatrix_wrp
  PUBLIC :: DistributedPairwiseMultiply_wrp
  PUBLIC :: DistributedGemm_wrp
  PUBLIC :: ScaleDistributedSparseMatrix_wrp
  PUBLIC :: DistributedSparseNorm_wrp
  PUBLIC :: Trace_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Matrix B = alpha*Matrix A + Matrix B (AXPY)
  SUBROUTINE IncrementDistributedSparseMatrix_wrp(ih_matA, ih_matB,&
       & alpha_in,threshold_in) &
       & bind(c,name="IncrementDistributedSparseMatrix_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matB(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: alpha_in
    REAL(NTREAL), INTENT(IN) :: threshold_in
    TYPE(DistributedSparseMatrix_wrp) :: h_matA
    TYPE(DistributedSparseMatrix_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    CALL IncrementDistributedSparseMatrix(h_matA%data, h_matB%data, alpha_in, &
         & threshold_in)
  END SUBROUTINE IncrementDistributedSparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> product = dot(matA,matB)
  FUNCTION DotDistributedSparseMatrix_wrp(ih_matA, ih_matB) RESULT(product)&
       & bind(c,name="DotDistributedSparseMatrix_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matB(SIZE_wrp)
    REAL(NTREAL) :: product
    TYPE(DistributedSparseMatrix_wrp) :: h_matA
    TYPE(DistributedSparseMatrix_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    product = DotDistributedSparseMatrix(h_matA%data, h_matB%data)
  END FUNCTION DotDistributedSparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Elementwise multiplication.
  SUBROUTINE DistributedPairwiseMultiply_wrp(ih_matA, ih_matB, ih_matC) &
       & bind(c,name="DistributedPairwiseMultiply_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_matB(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matC(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_matA
    TYPE(DistributedSparseMatrix_wrp) :: h_matB
    TYPE(DistributedSparseMatrix_wrp) :: h_matC

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    h_matC = TRANSFER(ih_matC,h_matC)

    CALL DistributedPairwiseMultiply(h_matA%data, h_matB%data, h_matC%data)
  END SUBROUTINE DistributedPairwiseMultiply_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Multiply two matrices together, and add to the third.
  SUBROUTINE DistributedGemm_wrp(ih_matA, ih_matB, ih_matC, alpha_in, &
       & beta_in, threshold_in, ih_memory_pool_in) &
       & bind(c,name="DistributedGemm_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_matB(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matC(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: alpha_in
    REAL(NTREAL), INTENT(IN) :: beta_in
    REAL(NTREAL), INTENT(IN) :: threshold_in
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_memory_pool_in(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_matA
    TYPE(DistributedSparseMatrix_wrp) :: h_matB
    TYPE(DistributedSparseMatrix_wrp) :: h_matC
    TYPE(DistributedMatrixMemoryPool_wrp) :: h_memory_pool_in

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    h_matC = TRANSFER(ih_matC,h_matC)
    h_memory_pool_in = TRANSFER(ih_memory_pool_in,h_memory_pool_in)

    CALL DistributedGemm(h_matA%data, h_matB%data, h_matC%data, &
         & alpha_in, beta_in, threshold_in, h_memory_pool_in%data)
  END SUBROUTINE DistributedGemm_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a distributed sparse matrix by a constant.
  SUBROUTINE ScaleDistributedSparseMatrix_wrp(ih_this, constant) &
       & bind(c,name="ScaleDistributedSparseMatrix_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: constant
    TYPE(DistributedSparseMatrix_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL ScaleDistributedSparseMatrix(h_this%data,constant)
  END SUBROUTINE ScaleDistributedSparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the norm of a distributed sparse matrix along the rows.
  FUNCTION DistributedSparseNorm_wrp(ih_this) &
       & bind(c,name="DistributedSparseNorm_wrp") RESULT(norm_value)
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    REAL(NTREAL) :: norm_value
    TYPE(DistributedSparseMatrix_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    norm_value = DistributedSparseNorm(h_this%data)
  END FUNCTION DistributedSparseNorm_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the trace of the matrix.
  FUNCTION Trace_wrp(ih_this) &
       & bind(c,name="Trace_wrp") RESULT(trace_value)
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    REAL(NTREAL) :: trace_value
    TYPE(DistributedSparseMatrix_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    trace_value = Trace(h_this%data)
  END FUNCTION Trace_wrp
END MODULE DistributedSparseMatrixAlgebraModule_wrp
