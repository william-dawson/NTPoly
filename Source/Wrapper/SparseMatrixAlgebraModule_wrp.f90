!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for wrapping the sparse algebra routines.
MODULE SparseMatrixAlgebraModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE MatrixMemoryPoolModule_wrp, ONLY : MatrixMemoryPool_wrp
  USE SparseMatrixAlgebraModule, ONLY : &
       & ScaleSparseMatrix, IncrementSparseMatrix, Gemm, &
       & DotSparseMatrix, PairwiseMultiplySparseMatrix
  USE SparseMatrixModule, ONLY : SparseMatrix_t
  USE SparseMatrixModule_wrp, ONLY: SparseMatrix_wrp
  USE TripletListModule_wrp, ONLY : TripletList_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int, c_char, c_bool
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ScaleSparseMatrix_wrp
  PUBLIC :: IncrementSparseMatrix_wrp
  PUBLIC :: PairwiseMultiplySparseMatrix_wrp
  PUBLIC :: DotSparseMatrix_wrp
  PUBLIC :: Gemm_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the scale a sparse matrix by a constant routine.
  !! @param[inout] ih_this handle to the matrix to scale.
  !! @param[in] constant scale factor.
  PURE SUBROUTINE ScaleSparseMatrix_wrp(ih_this, constant) &
       & bind(c,name="ScaleSparseMatrix_wrp")
    INTEGER(kind=c_int), INTENT(inout) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(in) :: constant
    TYPE(SparseMatrix_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL ScaleSparseMatrix(h_this%data,constant)
  END SUBROUTINE ScaleSparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap matrix incrementing function.
  !! @param[in] ih_matA handle to Matrix A.
  !! @param[in,out] ih_matB handle to Matrix B.
  !! @param[in] alpha_in multiplier.
  !! @param[in] threshold_in for flushing values to zero.
  PURE SUBROUTINE IncrementSparseMatrix_wrp(ih_matA, ih_matB, alpha_in, &
       & threshold_in) bind(c,name="IncrementSparseMatrix_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_matB(SIZE_wrp)
    REAL(NTREAL), INTENT(in) :: alpha_in
    REAL(NTREAL), INTENT(in) :: threshold_in
    TYPE(SparseMatrix_wrp) :: h_matA
    TYPE(SparseMatrix_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    CALL IncrementSparseMatrix(h_matA%data, h_matB%data, alpha_in, threshold_in)
  END SUBROUTINE IncrementSparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap matrix dot product function.
  !! @param[in] ih_matA handle to Matrix A.
  !! @param[in,out] ih_matB handle to Matrix B.
  PURE FUNCTION DotSparseMatrix_wrp(ih_matA, ih_matB) RESULT(product) &
       & bind(c,name="DotSparseMatrix_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: ih_matB(SIZE_wrp)
    REAL(NTREAL) :: product
    TYPE(SparseMatrix_wrp) :: h_matA
    TYPE(SparseMatrix_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    product = DotSparseMatrix(h_matA%data, h_matB%data)
  END FUNCTION DotSparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap pairwise matrix multiplication function.
  !! @param[in] ih_matA handle to Matrix A.
  !! @param[in] ih_matB handle to Matrix B.
  !! @param[out] ih_matC = A pairwise B
  SUBROUTINE PairwiseMultiplySparseMatrix_wrp(ih_matA, ih_matB, ih_matC) &
       & bind(c,name="PairwiseMultiplySparseMatrix_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: ih_matB(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_matC(SIZE_wrp)
    TYPE(SparseMatrix_wrp) :: h_matA
    TYPE(SparseMatrix_wrp) :: h_matB
    TYPE(SparseMatrix_wrp) :: h_matC

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    h_matC = TRANSFER(ih_matC,h_matC)

    CALL PairwiseMultiplySparseMatrix(h_matA%data, h_matB%data, h_matC%data)
  END SUBROUTINE PairwiseMultiplySparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap matrix multiplication function.
  !! @param[in] ih_matA handle to Matrix A.
  !! @param[in] ih_matB handle to Matrix B.
  !! @param[out] ih_matC = alpha*matA*op( matB ) + beta*matC.
  !! @param[in] IsATransposed true if A is already transposed.
  !! @param[in] IsBTransposed true if B is already transposed.
  !! @param[in] alpha scales the multiplication.
  !! @param[in] beta scales matrix we sum on to.
  !! @param[in] threshold for flushing values to zero.
  !! @param[inout] ih_blocked_memory_pool handle to memory pool.
  SUBROUTINE Gemm_wrp(ih_matA, ih_matB, ih_matC, IsATransposed, &
       & IsBTransposed, alpha, beta, threshold, ih_blocked_memory_pool) &
       & bind(c,name="Gemm_wrp")
    INTEGER(kind=c_int), INTENT(in) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(in) :: ih_matB(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(inout) :: ih_matC(SIZE_wrp)
    LOGICAL(kind=c_bool), INTENT(in) :: IsATransposed
    LOGICAL(kind=c_bool), INTENT(in) :: IsBTransposed
    REAL(NTREAL), INTENT(in) :: alpha
    REAL(NTREAL), INTENT(in) :: beta
    REAL(NTREAL), INTENT(in) :: threshold
    INTEGER(kind=c_int), INTENT(inout) :: ih_blocked_memory_pool(SIZE_wrp)
    TYPE(SparseMatrix_wrp) :: h_matA
    TYPE(SparseMatrix_wrp) :: h_matB
    TYPE(SparseMatrix_wrp) :: h_matC
    TYPE(MatrixMemoryPool_wrp) :: h_blocked_memory_pool

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    h_matC = TRANSFER(ih_matC,h_matC)
    h_blocked_memory_pool = TRANSFER(ih_blocked_memory_pool, &
         & h_blocked_memory_pool)

    CALL Gemm(h_matA%data, h_matB%data, h_matC%data, &
         & LOGICAL(IsATransposed), LOGICAL(IsBTransposed), alpha, &
         & beta, threshold, h_blocked_memory_pool%data)
  END SUBROUTINE Gemm_wrp
END MODULE SparseMatrixAlgebraModule_wrp
