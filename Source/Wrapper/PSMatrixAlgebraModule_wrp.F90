!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for wrapping a Distributed Sparse Matrix.
MODULE PSMatrixAlgebraModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE PMatrixMemoryPoolModule_wrp, ONLY : MatrixMemoryPool_p_wrp
  USE PSMatrixAlgebraModule
  USE PSMatrixModule_wrp, ONLY : Matrix_ps_wrp
  USE PermutationModule_wrp, ONLY : Permutation_wrp
  USE TripletListModule_wrp, ONLY : TripletList_r_wrp, TripletList_c_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int, c_char, c_bool
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: IncrementMatrix_ps_wrp
  PUBLIC :: DotMatrix_psr_wrp
  PUBLIC :: DotMatrix_psc_wrp
  PUBLIC :: MatrixPairwiseMultiply_ps_wrp
  PUBLIC :: MatrixMultiply_ps_wrp
  PUBLIC :: ScaleMatrix_ps_wrp
  PUBLIC :: MatrixNorm_ps_wrp
  PUBLIC :: MatrixTrace_ps_wrp
  PUBLIC :: MatrixDiagonalScale_psr_wrp
  PUBLIC :: MatrixDiagonalScale_psc_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Matrix B = alpha*Matrix A + Matrix B (AXPY)
  SUBROUTINE IncrementMatrix_ps_wrp(ih_matA, ih_matB, alpha_in,threshold_in) &
       & BIND(c,name="IncrementMatrix_ps_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matB(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: alpha_in
    REAL(NTREAL), INTENT(IN) :: threshold_in
    TYPE(Matrix_ps_wrp) :: h_matA
    TYPE(Matrix_ps_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    CALL IncrementMatrix(h_matA%DATA, h_matB%DATA, alpha_in, threshold_in)
  END SUBROUTINE IncrementMatrix_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> product = dot(matA,matB)
  SUBROUTINE DotMatrix_psr_wrp(ih_matA, ih_matB, product) &
       & BIND(c,name="DotMatrix_psr_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_matB(SIZE_wrp)
    REAL(NTREAL), INTENT(OUT) :: product
    TYPE(Matrix_ps_wrp) :: h_matA
    TYPE(Matrix_ps_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    CALL DotMatrix(h_matA%DATA, h_matB%DATA, product)
  END SUBROUTINE DotMatrix_psr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> product = dot(matA,matB)
  SUBROUTINE DotMatrix_psc_wrp(ih_matA, ih_matB, product_real, product_imag) &
       & BIND(c,name="DotMatrix_psc_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_matB(SIZE_wrp)
    REAL(NTREAL), INTENT(OUT) :: product_real
    REAL(NTREAL), INTENT(OUT) :: product_imag
    COMPLEX(NTREAL) :: product
    TYPE(Matrix_ps_wrp) :: h_matA
    TYPE(Matrix_ps_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    CALL DotMatrix(h_matA%DATA, h_matB%DATA, product)

    product_real = REAL(product)
    product_imag = AIMAG(product)
  END SUBROUTINE DotMatrix_psc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Elementwise multiplication.
  SUBROUTINE MatrixPairwiseMultiply_ps_wrp(ih_matA, ih_matB, ih_matC) &
       & BIND(c,name="MatrixPairwiseMultiply_ps_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_matB(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matC(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_matA
    TYPE(Matrix_ps_wrp) :: h_matB
    TYPE(Matrix_ps_wrp) :: h_matC

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    h_matC = TRANSFER(ih_matC,h_matC)

    CALL PairwiseMultiplyMatrix(h_matA%DATA, h_matB%DATA, h_matC%DATA)
  END SUBROUTINE MatrixPairwiseMultiply_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Multiply two matrices together, and add to the third.
  SUBROUTINE MatrixMultiply_ps_wrp(ih_matA, ih_matB, ih_matC, alpha_in, &
       & beta_in, threshold_in, ih_memory_pool_in) &
       & BIND(c,name="MatrixMultiply_ps_wrp")
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

    CALL MatrixMultiply(h_matA%DATA, h_matB%DATA, h_matC%DATA, &
         & alpha_in, beta_in, threshold_in, h_memory_pool_in%DATA)
  END SUBROUTINE MatrixMultiply_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a distributed sparse matrix by a constant.
  SUBROUTINE ScaleMatrix_ps_wrp(ih_this, constant) &
       & BIND(c,name="ScaleMatrix_ps_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: constant
    TYPE(Matrix_ps_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL ScaleMatrix(h_this%DATA,constant)
  END SUBROUTINE ScaleMatrix_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the norm of a distributed sparse matrix along the rows.
  FUNCTION MatrixNorm_ps_wrp(ih_this) BIND(c,name="MatrixNorm_ps_wrp") &
       & RESULT(norm_value)
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    REAL(NTREAL) :: norm_value
    TYPE(Matrix_ps_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    norm_value = MatrixNorm(h_this%DATA)
  END FUNCTION MatrixNorm_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the trace of the matrix.
  SUBROUTINE MatrixTrace_ps_wrp(ih_this, trace_value) &
       & BIND(c,name="MatrixTrace_ps_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(OUT) :: trace_value
    TYPE(Matrix_ps_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL MatrixTrace(h_this%DATA, trace_value)
  END SUBROUTINE MatrixTrace_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Scale a matrix using a diagonal matrix (triplet list form).
  SUBROUTINE MatrixDiagonalScale_psr_wrp(ih_mat, ih_tlist) &
       & BIND(c,name="MatrixDiagonalScale_psr_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_mat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_tlist(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_mat
    TYPE(TripletList_r_wrp) :: h_tlist

    h_mat = TRANSFER(ih_mat, h_mat)
    h_tlist = TRANSFER(ih_tlist, h_tlist)

    CALL MatrixDiagonalScale(h_mat%DATA, h_tlist%DATA)
  END SUBROUTINE MatrixDiagonalScale_psr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Scale a matrix using a diagonal matrix (triplet list form).
  SUBROUTINE MatrixDiagonalScale_psc_wrp(ih_mat, ih_tlist) &
       & BIND(c,name="MatrixDiagonalScale_psc_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_mat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_tlist(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_mat
    TYPE(TripletList_c_wrp) :: h_tlist

    h_mat = TRANSFER(ih_mat, h_mat)
    h_tlist = TRANSFER(ih_tlist, h_tlist)

    CALL MatrixDiagonalScale(h_mat%DATA, h_tlist%DATA)
  END SUBROUTINE MatrixDiagonalScale_psc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE PSMatrixAlgebraModule_wrp
