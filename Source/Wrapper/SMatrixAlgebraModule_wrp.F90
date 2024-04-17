!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for wrapping the sparse algebra routines.
MODULE SMatrixAlgebraModule_wrp
  USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
  USE MatrixMemoryPoolModule_wrp, ONLY : MatrixMemoryPool_lr_wrp, &
       & MatrixMemoryPool_lc_wrp
  USE SMatrixAlgebraModule
  USE SMatrixModule, ONLY : Matrix_lsr, Matrix_lsc
  USE SMatrixModule_wrp, ONLY: Matrix_lsr_wrp, Matrix_lsc_wrp
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
  PUBLIC :: DiagonalScale_lsr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ScaleMatrix_lsc_wrp
  PUBLIC :: IncrementMatrix_lsc_wrp
  PUBLIC :: DotMatrix_lsc_wrp
  PUBLIC :: PairwiseMultiplyMatrix_lsc_wrp
  PUBLIC :: MatrixMultiply_lsc_wrp
  PUBLIC :: DiagonalScale_lsc_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the scale a sparse matrix by a constant routine.
  SUBROUTINE ScaleMatrix_lsr_wrp(ih_this, constant) &
       & BIND(c,name="ScaleMatrix_lsr_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: constant
    TYPE(Matrix_lsr_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL ScaleMatrix(h_this%DATA,constant)
  END SUBROUTINE ScaleMatrix_lsr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap matrix incrementing function.
  SUBROUTINE IncrementMatrix_lsr_wrp(ih_matA, ih_matB, alpha_in, &
       & threshold_in) BIND(c,name="IncrementMatrix_lsr_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matB(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: alpha_in
    REAL(NTREAL), INTENT(IN) :: threshold_in
    TYPE(Matrix_lsr_wrp) :: h_matA
    TYPE(Matrix_lsr_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    CALL IncrementMatrix(h_matA%DATA, h_matB%DATA, alpha_in, threshold_in)
  END SUBROUTINE IncrementMatrix_lsr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap matrix dot product function.
  SUBROUTINE DotMatrix_lsr_wrp(ih_matA, ih_matB, product) &
       & BIND(c,name="DotMatrix_lsr_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_matB(SIZE_wrp)
    REAL(NTREAL), INTENT(OUT) :: product
    TYPE(Matrix_lsr_wrp) :: h_matA
    TYPE(Matrix_lsr_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    CALL DotMatrix(h_matA%DATA, h_matB%DATA, product)
  END SUBROUTINE DotMatrix_lsr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap pairwise matrix multiplication function.
  SUBROUTINE PairwiseMultiplyMatrix_lsr_wrp(ih_matA, ih_matB, ih_matC) &
       & BIND(c,name="PairwiseMultiplyMatrix_lsr_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_matB(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matC(SIZE_wrp)
    TYPE(Matrix_lsr_wrp) :: h_matA
    TYPE(Matrix_lsr_wrp) :: h_matB
    TYPE(Matrix_lsr_wrp) :: h_matC

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    h_matC = TRANSFER(ih_matC,h_matC)

    CALL PairwiseMultiplyMatrix(h_matA%DATA, h_matB%DATA, h_matC%DATA)
  END SUBROUTINE PairwiseMultiplyMatrix_lsr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap matrix multiplication function.
  SUBROUTINE MatrixMultiply_lsr_wrp(ih_matA, ih_matB, ih_matC, IsATransposed, &
       & IsBTransposed, alpha, beta, threshold, ih_blocked_memory_pool) &
       & BIND(c,name="MatrixMultiply_lsr_wrp")
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

    CALL MatrixMultiply(h_matA%DATA, h_matB%DATA, h_matC%DATA, &
         & LOGICAL(IsATransposed), LOGICAL(IsBTransposed), alpha, &
         & beta, threshold, h_blocked_memory_pool%DATA)
  END SUBROUTINE MatrixMultiply_lsr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Scale a matrix using a diagonal matrix (triplet list form).
  SUBROUTINE DiagonalScale_lsr_wrp(ih_mat, ih_tlist) &
       & BIND(c,name="DiagonalScale_lsr_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_mat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_tlist(SIZE_wrp)
    TYPE(Matrix_lsr_wrp) :: h_mat
    TYPE(TripletList_r_wrp) :: h_tlist

    h_mat = TRANSFER(ih_mat, h_mat)
    h_tlist = TRANSFER(ih_tlist, h_tlist)

    CALL DiagonalScale(h_mat%DATA, h_tlist%DATA)
  END SUBROUTINE DiagonalScale_lsr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the scale a sparse matrix by a constant routine.
  SUBROUTINE ScaleMatrix_lsc_wrp(ih_this, constant) &
       & BIND(c,name="ScaleMatrix_lsc_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: constant
    TYPE(Matrix_lsc_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL ScaleMatrix(h_this%DATA,constant)
  END SUBROUTINE ScaleMatrix_lsc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap matrix incrementing function.
  SUBROUTINE IncrementMatrix_lsc_wrp(ih_matA, ih_matB, alpha_in, &
       & threshold_in) BIND(c,name="IncrementMatrix_lsc_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matB(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: alpha_in
    REAL(NTREAL), INTENT(IN) :: threshold_in
    TYPE(Matrix_lsc_wrp) :: h_matA
    TYPE(Matrix_lsc_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    CALL IncrementMatrix(h_matA%DATA, h_matB%DATA, alpha_in, threshold_in)
  END SUBROUTINE IncrementMatrix_lsc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap matrix dot product function.
  SUBROUTINE DotMatrix_lsc_wrp(ih_matA, ih_matB, product_real, &
       & product_imag) BIND(c,name="DotMatrix_lsc_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_matB(SIZE_wrp)
    REAL(NTREAL), INTENT(OUT) :: product_real
    REAL(NTREAL), INTENT(OUT) :: product_imag
    COMPLEX(NTCOMPLEX) :: product
    TYPE(Matrix_lsc_wrp) :: h_matA
    TYPE(Matrix_lsc_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    CALL DotMatrix(h_matA%DATA, h_matB%DATA, product)

    product_real = REAL(product)
    product_imag = AIMAG(product)
  END SUBROUTINE DotMatrix_lsc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap pairwise matrix multiplication function.
  SUBROUTINE PairwiseMultiplyMatrix_lsc_wrp(ih_matA, ih_matB, ih_matC) &
       & BIND(c,name="PairwiseMultiplyMatrix_lsc_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_matB(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matC(SIZE_wrp)
    TYPE(Matrix_lsc_wrp) :: h_matA
    TYPE(Matrix_lsc_wrp) :: h_matB
    TYPE(Matrix_lsc_wrp) :: h_matC

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    h_matC = TRANSFER(ih_matC,h_matC)

    CALL PairwiseMultiplyMatrix(h_matA%DATA, h_matB%DATA, h_matC%DATA)
  END SUBROUTINE PairwiseMultiplyMatrix_lsc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap matrix multiplication function.
  SUBROUTINE MatrixMultiply_lsc_wrp(ih_matA, ih_matB, ih_matC, IsATransposed, &
       & IsBTransposed, alpha, beta, threshold, ih_blocked_memory_pool) &
       & BIND(c,name="MatrixMultiply_lsc_wrp")
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

    CALL MatrixMultiply(h_matA%DATA, h_matB%DATA, h_matC%DATA, &
         & LOGICAL(IsATransposed), LOGICAL(IsBTransposed), alpha, &
         & beta, threshold, h_blocked_memory_pool%DATA)
  END SUBROUTINE MatrixMultiply_lsc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Scale a matrix using a diagonal matrix (triplet list form).
  SUBROUTINE DiagonalScale_lsc_wrp(ih_mat, ih_tlist) &
       & BIND(c,name="DiagonalScale_lsc_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_mat(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_tlist(SIZE_wrp)
    TYPE(Matrix_lsc_wrp) :: h_mat
    TYPE(TripletList_c_wrp) :: h_tlist

    h_mat = TRANSFER(ih_mat, h_mat)
    h_tlist = TRANSFER(ih_tlist, h_tlist)

    CALL DiagonalScale(h_mat%DATA, h_tlist%DATA)
  END SUBROUTINE DiagonalScale_lsc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE SMatrixAlgebraModule_wrp
