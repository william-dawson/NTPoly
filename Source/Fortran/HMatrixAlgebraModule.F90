!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for hierarchy of blocked matrix types.
MODULE HMatrixAlgebraModule
  USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
  USE HMatrixModule, ONLY : Matrix_hsr, Matrix_hsc, GetMatrixBlockColumns, &
      & GetMatrixBlockRows, GetMatrixColumns, GetMatrixRows, IsBaseCase, &
      & ComposeMatrix, DestructMatrix, ConstructEmptyMatrix, IsAllocated
  USE SMatrixModule, ONLY : Matrix_lsr, Matrix_lsc
  USE SMatrixAlgebraModule, ONLY : ScaleMatrix, IncrementMatrix, DotMatrix, &
      & PairwiseMultiplyMatrix, MatrixMultiply, MatrixNorm, MatrixColumnNorm, &
      & MatrixGrandSum
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ScaleMatrix
  PUBLIC :: IncrementMatrix
  PUBLIC :: DotMatrix
  PUBLIC :: PairwiseMultiplyMatrix
  PUBLIC :: MatrixNorm
  PUBLIC :: MatrixColumnNorm
  PUBLIC :: MatrixGrandSum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE ScaleMatrix
     MODULE PROCEDURE ScaleMatrix_hsr
     MODULE PROCEDURE ScaleMatrix_hsc
     MODULE PROCEDURE ScaleMatrix_hsc_c
  END INTERFACE
  INTERFACE IncrementMatrix
     MODULE PROCEDURE IncrementMatrix_hsr
     MODULE PROCEDURE IncrementMatrix_hsc
  END INTERFACE
  INTERFACE DotMatrix
     MODULE PROCEDURE DotMatrix_hsr
     MODULE PROCEDURE DotMatrix_hsc
  END INTERFACE
  INTERFACE PairwiseMultiplyMatrix
     MODULE PROCEDURE PairwiseMultiplyMatrix_hsr
     MODULE PROCEDURE PairwiseMultiplyMatrix_hsc
  END INTERFACE
  INTERFACE MatrixNorm
     MODULE PROCEDURE MatrixNorm_hsr
     MODULE PROCEDURE MatrixNorm_hsc
  END INTERFACE
  INTERFACE MatrixColumnNorm
     MODULE PROCEDURE MatrixColumnNorm_hsr
     MODULE PROCEDURE MatrixColumnNorm_hsc
  END INTERFACE
  INTERFACE MatrixGrandSum
     MODULE PROCEDURE MatrixGrandSum_hsr
     MODULE PROCEDURE MatrixGrandSum_hsc
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a matrix by a constant.
  PURE SUBROUTINE ScaleMatrix_hsr(this, constant)
    !> The matrix to scale.
    TYPE(Matrix_hsr), INTENT(INOUT) :: this
    !> Constant scale factor.
    REAL(NTREAL), INTENT(IN) :: constant

    INCLUDE "block_algebra_includes/ScaleMatrix.f90"
  END SUBROUTINE ScaleMatrix_hsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a matrix by a constant.
  PURE SUBROUTINE ScaleMatrix_hsc(this, constant)
    !> The matrix to scale.
    TYPE(Matrix_hsc), INTENT(INOUT) :: this
    !> Constant scale factor.
    REAL(NTREAL), INTENT(IN) :: constant

    INCLUDE "block_algebra_includes/ScaleMatrix.f90"
  END SUBROUTINE ScaleMatrix_hsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a matrix by a constant.
  PURE SUBROUTINE ScaleMatrix_hsc_c(this, constant)
    !> The matrix to scale.
    TYPE(Matrix_hsc), INTENT(INOUT) :: this
    !> Constant scale factor.
    COMPLEX(NTCOMPLEX), INTENT(IN) :: constant

    INCLUDE "block_algebra_includes/ScaleMatrix.f90"
  END SUBROUTINE ScaleMatrix_hsc_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Matrix B = alpha*Matrix A + Matrix B (AXPY).
  PURE SUBROUTINE IncrementMatrix_hsr(matA, matB, alpha_in, threshold_in)
    !> Matrix A.
    TYPE(Matrix_hsr), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_hsr), INTENT(INOUT) :: matB
    !> Multiplier (default=1.0).
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    !> For flushing values to zero (default=0).
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    !! Local Variables
    TYPE(Matrix_hsr) :: matC

    INCLUDE "block_algebra_includes/IncrementMatrix.f90"
  END SUBROUTINE IncrementMatrix_hsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Matrix B = alpha*Matrix A + Matrix B (AXPY).
  PURE SUBROUTINE IncrementMatrix_hsc(matA, matB, alpha_in, threshold_in)
    !> Matrix A.
    TYPE(Matrix_hsc), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_hsc), INTENT(INOUT) :: matB
    !> Multiplier (default=1.0).
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    !> For flushing values to zero (default=0).
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    !! Local Variables
    TYPE(Matrix_hsc) :: matC

    INCLUDE "block_algebra_includes/IncrementMatrix.f90"
  END SUBROUTINE IncrementMatrix_hsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Product = sum(MatA[ij]*MatB[ij])
  PURE SUBROUTINE DotMatrix_hsr(matA, matB, product)
    !> Matrix A.
    TYPE(Matrix_hsr), INTENT(IN) :: matA
    !> Matrix B.
    TYPE(Matrix_hsr), INTENT(IN) :: matB
    !> Dot product.
    REAL(NTREAL), INTENT(OUT) :: product
    !! Local Variables
    REAL(NTREAL) :: temp

    INCLUDE "block_algebra_includes/DotMatrix.f90"
  END SUBROUTINE DotMatrix_hsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Product = sum(MatA^H[ij]*MatB[ij])
  PURE SUBROUTINE DotMatrix_hsc(matA, matB, product)
    !> Matrix A.
    TYPE(Matrix_hsc), INTENT(IN) :: matA
    !> Matrix B.
    TYPE(Matrix_hsc), INTENT(IN) :: matB
    !> Dot product.
    COMPLEX(NTCOMPLEX), INTENT(OUT) :: product
    !! Local Variables
    COMPLEX(NTCOMPLEX) :: temp

    INCLUDE "block_algebra_includes/DotMatrix.f90"
  END SUBROUTINE DotMatrix_hsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Pairwise Multiply two matrices.
  PURE SUBROUTINE PairwiseMultiplyMatrix_hsr(matA, matB, matC)
    !> Matrix A.
    TYPE(Matrix_hsr), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_hsr), INTENT(IN) :: matB
    !> matC = MatA mult MatB.
    TYPE(Matrix_hsr), INTENT(INOUT) :: matC

    INCLUDE "block_algebra_includes/PairwiseMultiplyMatrix.f90"
  END SUBROUTINE PairwiseMultiplyMatrix_hsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Pairwise Multiply two matrices.
  PURE SUBROUTINE PairwiseMultiplyMatrix_hsc(matA, matB, matC)
    !> Matrix A.
    TYPE(Matrix_hsc), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_hsc), INTENT(IN) :: matB
    !> matC = MatA mult MatB.
    TYPE(Matrix_hsc), INTENT(INOUT) :: matC

    INCLUDE "block_algebra_includes/PairwiseMultiplyMatrix.f90"
  END SUBROUTINE PairwiseMultiplyMatrix_hsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the 1 norm of a matrix.
  PURE FUNCTION MatrixNorm_hsr(this) RESULT(norm)
    !> The matrix to compute the norm of.
    TYPE(Matrix_hsr), INTENT(IN) :: this
    !> The norm of the matrix.
    REAL(NTREAL) :: norm

    INCLUDE "block_algebra_includes/MatrixNorm.f90"
  END FUNCTION MatrixNorm_hsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the 1 norm of a matrix.
  PURE FUNCTION MatrixNorm_hsc(this) RESULT(norm)
    !> The matrix to compute the norm of.
    TYPE(Matrix_hsc), INTENT(IN) :: this
    !> The norm of the matrix.
    REAL(NTREAL) :: norm

    INCLUDE "block_algebra_includes/MatrixNorm.f90"
  END FUNCTION MatrixNorm_hsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sum the elements of a matrix
  PURE SUBROUTINE MatrixGrandSum_hsr(this, sum_value)
    !> The matrix to sum
    TYPE(Matrix_hsr), INTENT(IN) :: this
    !> The sum of the matrix elements
    REAL(NTREAL), INTENT(OUT) :: sum_value
    !! Local Variables
    REAL(NTREAL) :: temp

    INCLUDE "block_algebra_includes/MatrixGrandSum.f90"
  END SUBROUTINE MatrixGrandSum_hsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sum the elements of a matrix
  PURE SUBROUTINE MatrixGrandSum_hsc(this, sum_value)
    !> The matrix to sum
    TYPE(Matrix_hsc), INTENT(IN) :: this
    !> The sum of the matrix elements
    COMPLEX(NTCOMPLEX), INTENT(OUT) :: sum_value
    !! Local Variables
    COMPLEX(NTCOMPLEX) :: temp

    INCLUDE "block_algebra_includes/MatrixGrandSum.f90"
  END SUBROUTINE MatrixGrandSum_hsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the norm of a sparse matrix along the columns.
  PURE SUBROUTINE MatrixColumnNorm_hsr(this, norm_per_column)
    !> The matrix to compute the norm of.
    TYPE(Matrix_hsr), INTENT(IN) :: this
    !> The norm value for each column in this matrix.
    REAL(NTREAL), DIMENSION(:), INTENT(OUT) :: norm_per_column
    !! Local Data
    TYPE(Matrix_lsr) :: merged
    REAL(NTREAL) :: temp_value

    INCLUDE "block_algebra_includes/MatrixColumnNorm.f90"
  END SUBROUTINE MatrixColumnNorm_hsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the norm of a sparse matrix along the columns.
  PURE SUBROUTINE MatrixColumnNorm_hsc(this, norm_per_column)
    !> The matrix to compute the norm of.
    TYPE(Matrix_hsc), INTENT(IN) :: this
    !> The norm value for each column in this matrix.
    REAL(NTREAL), DIMENSION(:), INTENT(OUT) :: norm_per_column
    !! Local Data
    TYPE(Matrix_lsc) :: merged
    COMPLEX(NTCOMPLEX)  :: temp_value

    INCLUDE "block_algebra_includes/MatrixColumnNorm.f90"
  END SUBROUTINE MatrixColumnNorm_hsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Multiply two matrices together, and add to the third.
  !> C := alpha*matA*op( matB ) + beta*matC
  SUBROUTINE MatrixMultiply_hsr(matA, matB, matC, IsATransposed_in, &
       & IsBTransposed_in, alpha_in, beta_in, threshold_in)
    !> Matrix A.
    TYPE(Matrix_hsr), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_hsr), INTENT(IN)  :: matB
    !> matC = alpha*matA*op( matB ) + beta*matC.
    TYPE(Matrix_hsr), INTENT(INOUT) :: matC
    !> True if blocks of A are already transposed.
    LOGICAL, OPTIONAL, INTENT(IN) :: IsATransposed_in
    !> True if blocks of B are already transposed.
    LOGICAL, OPTIONAL, INTENT(IN) :: IsBTransposed_in
    !> Scales the multiplication.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    !> Scales matrix we sum on to.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: beta_in
    !> For flushing values to zero. Default value is 0.0.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    !! Intermediate Data
    TYPE(Matrix_hsr) :: matAB
    LOGICAL :: IsATransposed, IsBTransposed
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: beta
    REAL(NTREAL) :: threshold

    INCLUDE "block_algebra_includes/MatrixMultiply.f90"
  END SUBROUTINE MatrixMultiply_hsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Multiply two matrices together, and add to the third.
  !> C := alpha*matA*op( matB ) + beta*matC
  SUBROUTINE MatrixMultiply_hsc(matA, matB, matC, IsATransposed_in, &
       & IsBTransposed_in, alpha_in, beta_in, threshold_in)
    !> Matrix A.
    TYPE(Matrix_hsc), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_hsc), INTENT(IN)  :: matB
    !> matC = alpha*matA*op( matB ) + beta*matC.
    TYPE(Matrix_hsc), INTENT(INOUT) :: matC
    !> True if blocks of A are already transposed.
    LOGICAL, OPTIONAL, INTENT(IN) :: IsATransposed_in
    !> True if blocks of B are already transposed.
    LOGICAL, OPTIONAL, INTENT(IN) :: IsBTransposed_in
    !> Scales the multiplication.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    !> Scales matrix we sum on to.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: beta_in
    !> For flushing values to zero. Default value is 0.0.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    !! Intermediate Data
    TYPE(Matrix_hsc) :: matAB
    LOGICAL :: IsATransposed, IsBTransposed
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: beta
    REAL(NTREAL) :: threshold

    INCLUDE "block_algebra_includes/MatrixMultiply.f90"
  END SUBROUTINE MatrixMultiply_hsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE HMatrixAlgebraModule
