!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for performing linear algebra using dense matrices.
MODULE MatrixDAlgebraModule
  USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
  USE MatrixDModule, ONLY : Matrix_ldr, Matrix_ldc
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ScaleMatrix
  PUBLIC :: IncrementMatrix
  PUBLIC :: DotMatrix
  PUBLIC :: PairwiseMultiplyMatrix
  PUBLIC :: MatrixMultiply
  PUBLIC :: MatrixNorm
  PUBLIC :: EigenDecomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE ScaleMatrix
     MODULE PROCEDURE ScaleMatrix_ldr
     MODULE PROCEDURE ScaleMatrix_ldc
  END INTERFACE
  INTERFACE IncrementMatrix
     MODULE PROCEDURE IncrementMatrix_ldr
     MODULE PROCEDURE IncrementMatrix_ldc
  END INTERFACE
  INTERFACE DotMatrix
     MODULE PROCEDURE DotMatrix_ldr
     MODULE PROCEDURE DotMatrix_ldc
  END INTERFACE
  INTERFACE PairwiseMultiplyMatrix
     MODULE PROCEDURE PairwiseMultiplyMatrix_ldr
     MODULE PROCEDURE PairwiseMultiplyMatrix_ldc
  END INTERFACE
  INTERFACE MatrixMultiply
     MODULE PROCEDURE MatrixMultiply_ldr
     MODULE PROCEDURE MatrixMultiply_ldc
  END INTERFACE
  INTERFACE MatrixNorm
     MODULE PROCEDURE MatrixNorm_ldr
     MODULE PROCEDURE MatrixNorm_ldc
  END INTERFACE
  INTERFACE EigenDecomposition
     MODULE PROCEDURE EigenDecomposition_ldr
     MODULE PROCEDURE EigenDecomposition_ldc
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a sparse matrix by a constant.
  !! @param[inout] matA Matrix A.
  !! @param[in] constant scale factor.
  PURE SUBROUTINE ScaleMatrix_ldr(matA,constant)
    !! Parameters
    TYPE(Matrix_ldr), INTENT(INOUT) :: matA
    REAL(NTREAL), INTENT(IN) :: constant

    matA%data = matA%data*constant
  END SUBROUTINE ScaleMatrix_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a sparse matrix by a constant.
  !! @param[inout] matA Matrix A.
  !! @param[in] constant scale factor.
  PURE SUBROUTINE ScaleMatrix_ldc(matA,constant)
    !! Parameters
    TYPE(Matrix_ldc), INTENT(INOUT) :: matA
    REAL(NTREAL), INTENT(IN) :: constant

    matA%data = matA%data*constant
  END SUBROUTINE ScaleMatrix_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Matrix B = alpha*Matrix A + Matrix B (AXPY).
  !! This will utilize the sparse vector addition routine.
  !! @param[in] matA Matrix A.
  !! @param[in,out] matB Matrix B.
  !! @param[in] alpha_in multiplier (optional, default=1.0)
  PURE SUBROUTINE IncrementMatrix_ldr(matA, matB, alpha_in)
    !! Parameters
    TYPE(Matrix_ldr), INTENT(IN)  :: matA
    TYPE(Matrix_ldr), INTENT(INOUT) :: matB
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in

    INCLUDE "dense_includes/IncrementMatrix.f90"
  END SUBROUTINE IncrementMatrix_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Matrix B = alpha*Matrix A + Matrix B (AXPY).
  !! This will utilize the sparse vector addition routine.
  !! @param[in] matA Matrix A.
  !! @param[in,out] matB Matrix B.
  !! @param[in] alpha_in multiplier (optional, default=1.0)
  PURE SUBROUTINE IncrementMatrix_ldc(matA, matB, alpha_in)
    !! Parameters
    TYPE(Matrix_ldc), INTENT(IN)  :: matA
    TYPE(Matrix_ldc), INTENT(INOUT) :: matB
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in

    INCLUDE "dense_includes/IncrementMatrix.f90"
  END SUBROUTINE IncrementMatrix_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Product = sum(MatA[ij]*MatB[ij])
  !! @param[in] matA Matrix A.
  !! @param[in] matB Matrix B.
  !! @result product
  PURE FUNCTION DotMatrix_ldr(matA, matB) RESULT(product)
    !! Parameters
    TYPE(Matrix_ldr), INTENT(IN) :: matA
    TYPE(Matrix_ldr), INTENT(IN) :: matB
    REAL(NTREAL) :: product
    !! Local Variables
    TYPE(Matrix_ldr) :: matC

    CALL PairwiseMultiplyMatrix(matA, matB, matC)
    product = SUM(matC%data)
    CALL matC%Destruct

  END FUNCTION DotMatrix_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Product = sum(MatA[ij]*MatB[ij])
  !! @param[in] matA Matrix A.
  !! @param[in] matB Matrix B.
  !! @result product
  PURE FUNCTION DotMatrix_ldc(matA, matB) RESULT(product)
    !! Parameters
    TYPE(Matrix_ldc), INTENT(IN) :: matA
    TYPE(Matrix_ldc), INTENT(IN) :: matB
    REAL(NTREAL) :: product
    !! Local Variables
    TYPE(Matrix_ldc) :: matC

    CALL PairwiseMultiplyMatrix(matA, matB, matC)
    product = SUM(matC%data)
    CALL matC%Destruct

  END FUNCTION DotMatrix_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Pairwise Multiply two matrices.
  !! This will utilize the sparse vector pairwise routine.
  !! @param[in] matA Matrix A.
  !! @param[in] matB Matrix B.
  !! @param[in,out] matC = MatA mult MatB.
  PURE SUBROUTINE PairwiseMultiplyMatrix_ldr(matA, matB, matC)
    !! Parameters
    TYPE(Matrix_ldr), INTENT(IN)  :: matA
    TYPE(Matrix_ldr), INTENT(IN) :: matB
    TYPE(Matrix_ldr), INTENT(INOUT) :: matC

    INCLUDE "dense_includes/PairwiseMultiplyMatrix.f90"
  END SUBROUTINE PairwiseMultiplyMatrix_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Pairwise Multiply two matrices.
  !! This will utilize the sparse vector pairwise routine.
  !! @param[in] matA Matrix A.
  !! @param[in] matB Matrix B.
  !! @param[in,out] matC = MatA mult MatB.
  PURE SUBROUTINE PairwiseMultiplyMatrix_ldc(matA, matB, matC)
    !! Parameters
    TYPE(Matrix_ldc), INTENT(IN)  :: matA
    TYPE(Matrix_ldc), INTENT(IN) :: matB
    TYPE(Matrix_ldc), INTENT(INOUT) :: matC

    INCLUDE "dense_includes/PairwiseMultiplyMatrix.f90"
  END SUBROUTINE PairwiseMultiplyMatrix_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for multiplying two dense matrices.
  !! @param[in] MatA the first matrix.
  !! @param[in] MatB the second matrix.
  !! @param[inout] MatC = MatA*MatB.
  SUBROUTINE MatrixMultiply_ldr(MatA,MatB,MatC)
    !! Parameters
    TYPE(Matrix_ldr), INTENT(IN) :: MatA
    TYPE(Matrix_ldr), INTENT(IN) :: MatB
    TYPE(Matrix_ldr), INTENT(INOUT) :: MatC
    !! Local variables
    CHARACTER, PARAMETER :: TRANSA = 'N'
    CHARACTER, PARAMETER :: TRANSB = 'N'
    INTEGER :: M
    INTEGER :: N
    INTEGER :: K
    DOUBLE PRECISION, PARAMETER :: ALPHA = 1.0
    INTEGER :: LDA
    INTEGER :: LDB
    DOUBLE PRECISION, PARAMETER :: BETA = 0.0
    INTEGER :: LDC

    CALL MatC%InitEmpty(MatA%rows, MatB%columns)

    !! Setup Lapack
    M = MatA%rows
    N = MatB%columns
    K = MatA%columns
    LDA = M
    LDB = K
    LDC = M

    CALL DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, MatA%data, LDA, MatB%data, &
         & LDB, BETA, MatC%data, LDC)

  END SUBROUTINE MatrixMultiply_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for multiplying two dense matrices.
  !! @param[in] MatA the first matrix.
  !! @param[in] MatB the second matrix.
  !! @param[inout] MatC = MatA*MatB.
  SUBROUTINE MatrixMultiply_ldc(MatA,MatB,MatC)
    !! Parameters
    TYPE(Matrix_ldc), INTENT(IN) :: MatA
    TYPE(Matrix_ldc), INTENT(IN) :: MatB
    TYPE(Matrix_ldc), INTENT(INOUT) :: MatC
    !! Local variables
    CHARACTER, PARAMETER :: TRANSA = 'N'
    CHARACTER, PARAMETER :: TRANSB = 'N'
    INTEGER :: M
    INTEGER :: N
    INTEGER :: K
    COMPLEX*16, PARAMETER :: ALPHA = 1.0
    INTEGER :: LDA
    INTEGER :: LDB
    COMPLEX*16, PARAMETER :: BETA = 0.0
    INTEGER :: LDC

    CALL MatC%InitEmpty(MatA%rows, MatB%columns)

    !! Setup Lapack
    M = MatA%rows
    N = MatB%columns
    K = MatA%columns
    LDA = M
    LDB = K
    LDC = M

    CALL ZGEMM(TRANSA, TRANSB, M, N, K, ALPHA, MatA%data, LDA, MatB%data, &
         & LDB, BETA, MatC%data, LDC)

  END SUBROUTINE MatrixMultiply_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the norm of a dense matrix.
  !! Computes the Frobenius norm.
  !! @param[in] this the matrix to compute the norm of.
  !! @result the norm of the matrix.
  FUNCTION MatrixNorm_ldr(this) RESULT(norm)
    !! Parameters
    TYPE(Matrix_ldr), INTENT(IN) :: this
    REAL(NTREAL) :: norm
    !! Local Variables
    CHARACTER, PARAMETER :: NORMC = 'F'
    DOUBLE PRECISION, DIMENSION(this%rows) :: WORK
    INTEGER :: M, N, LDA
    !! Externel
    DOUBLE PRECISION, EXTERNAL :: dlange

    M = this%rows
    N = this%columns
    LDA = this%rows
    norm = DLANGE(NORMC, M, N, this%data, LDA, WORK)
  END FUNCTION MatrixNorm_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the norm of a dense matrix.
  !! Computes the Frobenius norm.
  !! @param[in] this the matrix to compute the norm of.
  !! @result the norm of the matrix.
  FUNCTION MatrixNorm_ldc(this) RESULT(norm)
    !! Parameters
    TYPE(Matrix_ldc), INTENT(IN) :: this
    REAL(NTREAL) :: norm
    !! Local Variables
    CHARACTER, PARAMETER :: NORMC = 'F'
    COMPLEX*16, DIMENSION(this%rows) :: WORK
    INTEGER :: M, N, LDA
    !! Externel
    DOUBLE PRECISION, EXTERNAL :: zlange

    M = this%rows
    N = this%columns
    LDA = this%rows
    norm = ZLANGE(NORMC, M, N, this%data, LDA, WORK)
  END FUNCTION MatrixNorm_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the eigenvectors of a dense matrix.
  !! Wraps a standard dense linear algebra routine.
  !! @param[in] MatA the matrix to decompose.
  !! @param[out] MatV the eigenvectors.
  !! @param[out] MatV the eigenvalues.
  SUBROUTINE EigenDecomposition_ldr(MatA, MatV, MatW)
    !! Parameters
    TYPE(Matrix_ldr), INTENT(IN) :: MatA
    TYPE(Matrix_ldr), INTENT(INOUT) :: MatV
    TYPE(Matrix_ldr), INTENT(INOUT), OPTIONAL :: MatW
    !! Local variables
    CHARACTER, PARAMETER :: job = 'V', uplo = 'U'
    INTEGER :: N, LDA
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: W
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WORK
    DOUBLE PRECISION :: WORKTEMP
    INTEGER :: LWORK
    INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
    INTEGER :: IWORKTEMP
    INTEGER :: LIWORK
    INTEGER :: INFO
    INTEGER :: II

    CALL MatV%InitEmpty(MatA%rows, MatA%columns)
    MatV%data = MatA%data

    N = SIZE(MatA%data,DIM=1)
    LDA = N

    !! Allocations
    ALLOCATE(W(N))

    !! Determine the scratch space size
    LWORK = -1
    CALL DSYEVD(JOB, UPLO, N, MatA%data, LDA, W, WORKTEMP, LWORK, IWORKTEMP, &
         & LIWORK, INFO)
    N = LDA
    LWORK = INT(WORKTEMP)
    ALLOCATE(WORK(LWORK))
    LIWORK = INT(IWORKTEMP)
    ALLOCATE(IWORK(LIWORK))

    !! Run Lapack For Real
    CALL DSYEVD(JOB, UPLO, N, MatV%data, LDA, W, WORK, LWORK, IWORK, LIWORK, &
         & INFO)

    !! Extract Eigenvalues
    IF (PRESENT(MatW)) THEN
       CALL MatW%InitEmpty(MatA%rows,1)
       DO II = 1, N
          MatW%data(II,1) = W(II)
       END DO
    END IF

    !! Cleanup
    DEALLOCATE(W)
    DEALLOCATE(Work)

  END SUBROUTINE EigenDecomposition_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the eigenvectors of a dense matrix.
  !! Wraps a standard dense linear algebra routine.
  !! @param[in] MatA the matrix to decompose.
  !! @param[out] MatV the eigenvectors.
  !! @param[out] MatV the eigenvalues.
  SUBROUTINE EigenDecomposition_ldc(MatA, MatV, MatW)
    !! Parameters
    TYPE(Matrix_ldc), INTENT(IN) :: MatA
    TYPE(Matrix_ldc), INTENT(INOUT) :: MatV
    TYPE(Matrix_ldc), INTENT(INOUT), OPTIONAL :: MatW
    !! Standard parameters
    CHARACTER, PARAMETER :: job = 'V', uplo = 'U'
    INTEGER :: N, LDA
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: W
    COMPLEX*16, DIMENSION(:), ALLOCATABLE :: WORK
    INTEGER :: LWORK
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RWORK
    INTEGER :: LRWORK
    INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
    INTEGER :: LIWORK
    INTEGER :: INFO
    !! Temp
    COMPLEX*16 :: WORKTEMP
    DOUBLE PRECISION :: RWORKTEMP
    INTEGER :: IWORKTEMP
    INTEGER :: II

    CALL MatV%InitEmpty(MatA%rows, MatA%columns)
    MatV%data = MatA%data

    N = SIZE(MatA%data,DIM=1)
    LDA = N

    !! Allocations
    ALLOCATE(W(N))

    !! Determine the scratch space size
    LWORK = -1
    CALL ZHEEVD(JOB, UPLO, N, MatA%data, LDA, W, WORKTEMP, LWORK, RWORKTEMP, &
         & LRWORK, IWORKTEMP, LIWORK, INFO)
    N = LDA
    LWORK = INT(WORKTEMP)
    ALLOCATE(WORK(LWORK))
    LRWORK = INT(RWORKTEMP)
    ALLOCATE(RWORK(LRWORK))
    LIWORK = INT(IWORKTEMP)
    ALLOCATE(IWORK(LIWORK))

    !! Run Lapack For Real
    CALL ZHEEVD(JOB, UPLO, N, MatV%data, LDA, W, WORK, LWORK, RWORK, LRWORK, &
         & IWORK, LIWORK, INFO)

    !! Extract Eigenvalues
    IF (PRESENT(MatW)) THEN
       CALL MatW%InitEmpty(MatA%rows, 1)
       DO II = 1, N
          MatW%data(II,1) = W(II)
       END DO
    END IF

    !! Cleanup
    DEALLOCATE(W)
    DEALLOCATE(Work)
    DEALLOCATE(RWork)

  END SUBROUTINE EigenDecomposition_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE MatrixDAlgebraModule
