!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Performing Distributed Sparse Matrix Algebra Operations.
MODULE PSMatrixAlgebraModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL, NTCOMPLEX, MPINTCOMPLEX
  USE GemmTasksModule
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & CheckMemoryPoolValidity, DestructMatrixMemoryPool
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, CopyMatrix, &
       & DestructMatrix, ConvertMatrixToComplex, ConjugateMatrix, &
       & MergeMatrixLocalBlocks
  USE MatrixReduceModule, ONLY : ReduceHelper_t, ReduceAndComposeMatrixSizes, &
       & ReduceAndComposeMatrixData, ReduceAndComposeMatrixCleanup, &
       & ReduceANdSumMatrixSizes, ReduceAndSumMatrixData, &
       & ReduceAndSumMatrixCleanup, TestReduceSizeRequest, &
       & TestReduceInnerRequest, TestReduceDataRequest
  USE SMatrixAlgebraModule, ONLY : MatrixMultiply, MatrixGrandSum, &
       & PairwiseMultiplyMatrix, IncrementMatrix, ScaleMatrix, &
       & MatrixColumnNorm
  USE SMatrixModule, ONLY : Matrix_lsr, Matrix_lsc, DestructMatrix, CopyMatrix,&
       & TransposeMatrix, ComposeMatrixColumns, MatrixToTripletList
  USE TimerModule, ONLY : StartTimer, StopTimer
  USE TripletListModule, ONLY : TripletList_r, TripletList_c
  USE NTMPIModule
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: MatrixSigma
  PUBLIC :: MatrixMultiply
  PUBLIC :: MatrixGrandSum
  PUBLIC :: PairwiseMultiplyMatrix
  PUBLIC :: MatrixNorm
  PUBLIC :: DotMatrix
  PUBLIC :: IncrementMatrix
  PUBLIC :: ScaleMatrix
  PUBLIC :: MatrixTrace
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE MatrixSigma
     MODULE PROCEDURE MatrixSigma_ps
  END INTERFACE
  INTERFACE MatrixMultiply
     MODULE PROCEDURE MatrixMultiply_ps
  END INTERFACE
  INTERFACE MatrixGrandSum
     MODULE PROCEDURE MatrixGrandSum_psr
     MODULE PROCEDURE MatrixGrandSum_psc
  END INTERFACE
  INTERFACE PairwiseMultiplyMatrix
     MODULE PROCEDURE PairwiseMultiplyMatrix_ps
  END INTERFACE
  INTERFACE MatrixNorm
     MODULE PROCEDURE MatrixNorm_ps
  END INTERFACE
  INTERFACE DotMatrix
     MODULE PROCEDURE DotMatrix_psr
     MODULE PROCEDURE DotMatrix_psc
  END INTERFACE
  INTERFACE IncrementMatrix
     MODULE PROCEDURE IncrementMatrix_ps
  END INTERFACE
  INTERFACE ScaleMatrix
     MODULE PROCEDURE ScaleMatrix_psr
     MODULE PROCEDURE ScaleMatrix_psc
  END INTERFACE
  INTERFACE MatrixTrace
     MODULE PROCEDURE MatrixTrace_psr
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute sigma for the inversion method.
  !> See \cite ozaki2001efficient for details.
  SUBROUTINE MatrixSigma_ps(this, sigma_value)
    !> The matrix to compute the sigma value of.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> Sigma
    REAL(NTREAL), INTENT(OUT) :: sigma_value
    !! Local Data
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: column_sigma_contribution
    !! Counters/Temporary
    INTEGER :: inner_counter, outer_counter
    TYPE(Matrix_lsr) :: merged_local_data_r
    TYPE(Matrix_lsc) :: merged_local_data_c
    INTEGER :: ierr

    IF (this%is_complex) THEN
#define LMAT merged_local_data_c
#include "distributed_algebra_includes/MatrixSigma.f90"
#undef LMAT
    ELSE
#define LMAT merged_local_data_r
#include "distributed_algebra_includes/MatrixSigma.f90"
#undef LMAT
    ENDIF
  END SUBROUTINE MatrixSigma_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Multiply two matrices together, and add to the third.
  !> C := alpha*matA*matB+ beta*matC
  SUBROUTINE MatrixMultiply_ps(matA, matB ,matC, alpha_in, beta_in, &
       & threshold_in, memory_pool_in)
    !> Matrix A.
    TYPE(Matrix_ps), INTENT(IN)        :: matA
    !> Matrix B.
    TYPE(Matrix_ps), INTENT(IN)        :: matB
    !> matC = alpha*matA*matB + beta*matC
    TYPE(Matrix_ps), INTENT(INOUT)     :: matC
    !> Scales the multiplication
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    !> Scales matrix we sum on to.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: beta_in
    !> For flushing values to zero.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    !> A memory pool for the calculation.
    TYPE(MatrixMemoryPool_p), OPTIONAL, INTENT(INOUT) :: memory_pool_in
    !! Local Versions of Optional Parameter
    TYPE(Matrix_ps) :: matAConverted
    TYPE(Matrix_ps) :: matBConverted
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: beta
    REAL(NTREAL) :: threshold
    TYPE(MatrixMemoryPool_p) :: memory_pool

    !! Handle the optional parameters
    IF (.NOT. PRESENT(alpha_in)) THEN
       alpha = 1.0d+0
    ELSE
       alpha = alpha_in
    END IF
    IF (.NOT. PRESENT(beta_in)) THEN
       beta = 0.0
    ELSE
       beta = beta_in
    END IF
    IF (.NOT. PRESENT(threshold_in)) THEN
       threshold = 0.0
    ELSE
       threshold = threshold_in
    END IF

    !! Setup Memory Pool
    IF (PRESENT(memory_pool_in)) THEN
       IF (matA%is_complex) THEN
          IF (.NOT. CheckMemoryPoolValidity(memory_pool_in, matA)) THEN
             CALL DestructMatrixMemoryPool(memory_pool_in)
             memory_pool_in = MatrixMemoryPool_p(matA)
          END IF
       ELSE
          IF (.NOT. CheckMemoryPoolValidity(memory_pool_in, matB)) THEN
             CALL DestructMatrixMemoryPool(memory_pool_in)
             memory_pool_in = MatrixMemoryPool_p(matB)
          END IF
       END IF
    ELSE
       IF (matA%is_complex) THEN
          memory_pool = MatrixMemoryPool_p(matA)
       ELSE
          memory_pool = MatrixMemoryPool_p(matB)
       END IF
    END IF


    !! Perform Upcasting
    IF (matB%is_complex .AND. .NOT. matA%is_complex) THEN
       CALL ConvertMatrixToComplex(matA, matAConverted)
       IF (PRESENT(memory_pool_in)) THEN
          CALL MatrixMultiply_ps_imp(matAConverted, matB, matC, alpha, beta, &
               & threshold, memory_pool_in)
       ELSE
          CALL MatrixMultiply_ps_imp(matAConverted, matB, matC, alpha, beta, &
               & threshold, memory_pool)
       END IF
    ELSE IF (matA%is_complex .AND. .NOT. matB%is_complex) THEN
       CALL ConvertMatrixToComplex(matB, matBConverted)
       IF (PRESENT(memory_pool_in)) THEN
          CALL MatrixMultiply_ps_imp(matA, matBConverted, matC, alpha, beta, &
               & threshold, memory_pool_in)
       ELSE
          CALL MatrixMultiply_ps_imp(matA, matBConverted, matC, alpha, beta, &
               & threshold, memory_pool)
       END IF
    ELSE
       IF (PRESENT(memory_pool_in)) THEN
          CALL MatrixMultiply_ps_imp(matA, matB, matC, alpha, beta, &
               & threshold, memory_pool_in)
       ELSE
          CALL MatrixMultiply_ps_imp(matA, matB, matC, alpha, beta, &
               & threshold, memory_pool)
       END IF
    END IF

    CALL DestructMatrixMemoryPool(memory_pool)
    CALL DestructMatrix(matAConverted)
    CALL DestructMatrix(matBConverted)

  END SUBROUTINE MatrixMultiply_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The actual implementation of matrix multiply is here. Takes the
  !> same parameters as the standard multiply, but nothing is optional.
  SUBROUTINE MatrixMultiply_ps_imp(matA, matB ,matC, alpha, beta, &
       & threshold, memory_pool)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN)    :: matA
    TYPE(Matrix_ps), INTENT(IN)    :: matB
    TYPE(Matrix_ps), INTENT(INOUT) :: matC
    REAL(NTREAL), INTENT(IN) :: alpha
    REAL(NTREAL), INTENT(IN) :: beta
    REAL(NTREAL), INTENT(IN) :: threshold
    TYPE(MatrixMemoryPool_p), INTENT(INOUT) :: memory_pool
    TYPE(Matrix_ps) :: matAB
    !! Temporary Matrices
    TYPE(Matrix_lsr), DIMENSION(:,:), ALLOCATABLE :: AdjacentABlocks_r
    TYPE(Matrix_lsr), DIMENSION(:), ALLOCATABLE :: LocalRowContribution_r
    TYPE(Matrix_lsr), DIMENSION(:), ALLOCATABLE :: GatheredRowContribution_r
    TYPE(Matrix_lsr), DIMENSION(:), ALLOCATABLE :: GatheredRowContributionT_r
    TYPE(Matrix_lsr), DIMENSION(:,:), ALLOCATABLE :: TransposedBBlocks_r
    TYPE(Matrix_lsr), DIMENSION(:), ALLOCATABLE :: LocalColumnContribution_r
    TYPE(Matrix_lsr), DIMENSION(:), ALLOCATABLE :: GatheredColumnContribution_r
    TYPE(Matrix_lsr), DIMENSION(:,:), ALLOCATABLE :: SliceContribution_r
    TYPE(Matrix_lsc), DIMENSION(:,:), ALLOCATABLE :: AdjacentABlocks_c
    TYPE(Matrix_lsc), DIMENSION(:), ALLOCATABLE :: LocalRowContribution_c
    TYPE(Matrix_lsc), DIMENSION(:), ALLOCATABLE :: GatheredRowContribution_c
    TYPE(Matrix_lsc), DIMENSION(:), ALLOCATABLE :: GatheredRowContributionT_c
    TYPE(Matrix_lsc), DIMENSION(:,:), ALLOCATABLE :: TransposedBBlocks_c
    TYPE(Matrix_lsc), DIMENSION(:), ALLOCATABLE :: LocalColumnContribution_c
    TYPE(Matrix_lsc), DIMENSION(:), ALLOCATABLE :: GatheredColumnContribution_c
    TYPE(Matrix_lsc), DIMENSION(:,:), ALLOCATABLE :: SliceContribution_c
    !! Communication Helpers
    TYPE(ReduceHelper_t), DIMENSION(:), ALLOCATABLE :: row_helper
    TYPE(ReduceHelper_t), DIMENSION(:), ALLOCATABLE :: column_helper
    TYPE(ReduceHelper_t), DIMENSION(:,:), ALLOCATABLE :: slice_helper
    !! For Iterating Over Local Blocks
    INTEGER :: II, II2
    INTEGER :: JJ, JJ2
    INTEGER :: duplicate_start_column, duplicate_offset_column
    INTEGER :: duplicate_start_row, duplicate_offset_row
    REAL(NTREAL) :: working_threshold
    !! Scheduling the A work
    INTEGER, DIMENSION(:), ALLOCATABLE :: ATasks
    INTEGER :: ATasks_completed
    !! Scheduling the B work
    INTEGER, DIMENSION(:), ALLOCATABLE :: BTasks
    INTEGER :: BTasks_completed
    !! Scheduling the AB work
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ABTasks
    INTEGER :: ABTasks_completed

    IF (matA%is_complex) THEN
#define AdjacentABlocks AdjacentABlocks_c
#define LocalRowContribution LocalRowContribution_c
#define GatheredRowContribution GatheredRowContribution_c
#define GatheredRowContributionT GatheredRowContributionT_c
#define TransposedBBlocks TransposedBBlocks_c
#define LocalColumnContribution LocalColumnContribution_c
#define GatheredColumnContribution GatheredColumnContribution_c
#define SliceContribution SliceContribution_c
#define LMAT local_data_c
#define MPGRID memory_pool%grid_c
#include "distributed_algebra_includes/MatrixMultiply.f90"
#undef AdjacentABlocks
#undef LocalRowContribution
#undef GatheredRowContribution
#undef GatheredRowContributionT
#undef TransposedBBlocks
#undef LocalColumnContribution
#undef GatheredColumnContribution
#undef SliceContribution
#undef LMAT
#undef MPGRID
    ELSE
#define AdjacentABlocks AdjacentABlocks_r
#define LocalRowContribution LocalRowContribution_r
#define GatheredRowContribution GatheredRowContribution_r
#define GatheredRowContributionT GatheredRowContributionT_r
#define TransposedBBlocks TransposedBBlocks_r
#define LocalColumnContribution LocalColumnContribution_r
#define GatheredColumnContribution GatheredColumnContribution_r
#define SliceContribution SliceContribution_r
#define LMAT local_data_r
#define MPGRID memory_pool%grid_r
#include "distributed_algebra_includes/MatrixMultiply.f90"
#undef AdjacentABlocks
#undef LocalRowContribution
#undef GatheredRowContribution
#undef GatheredRowContributionT
#undef TransposedBBlocks
#undef LocalColumnContribution
#undef GatheredColumnContribution
#undef SliceContribution
#undef LMAT
#undef MPGRID
    END IF
  END SUBROUTINE MatrixMultiply_ps_imp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sum up the elements in a matrix into a single value.
  SUBROUTINE MatrixGrandSum_psr(this, sum)
    !> The matrix to compute.
    TYPE(Matrix_ps), INTENT(IN)  :: this
    !> The sum of all elements.
    REAL(NTREAL), INTENT(OUT) :: sum
    !! Local Data
    INTEGER :: II, JJ
    REAL(NTREAL) :: temp_r
    COMPLEX(NTCOMPLEX) :: temp_c
    INTEGER :: ierr

#define MPIDATATYPE MPINTREAL
    IF (this%is_complex) THEN
#define TEMP temp_c
#define LMAT local_data_c
#include "distributed_algebra_includes/MatrixGrandSum.f90"
#undef LMAT
#undef TEMP
    ELSE
#define TEMP temp_r
#define LMAT local_data_r
#include "distributed_algebra_includes/MatrixGrandSum.f90"
#undef LMAT
#undef TEMP
    END IF
#undef MPIDATATYPE
  END SUBROUTINE MatrixGrandSum_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sum up the elements in a matrix into a single value.
  SUBROUTINE MatrixGrandSum_psc(this, sum)
    !! Parameters
    !> The matrix to compute.
    TYPE(Matrix_ps), INTENT(IN)  :: this
    !> The sum of all elements.
    COMPLEX(NTCOMPLEX), INTENT(OUT) :: sum
    !! Local Data
    INTEGER :: II, JJ
    REAL(NTREAL) :: temp_r
    COMPLEX(NTCOMPLEX) :: temp_c
    INTEGER :: ierr

#define MPIDATATYPE MPINTCOMPLEX
    IF (this%is_complex) THEN
#define TEMP temp_c
#define LMAT local_data_c
#include "distributed_algebra_includes/MatrixGrandSum.f90"
#undef LMAT
#undef TEMP
    ELSE
#define TEMP temp_r
#define LMAT local_data_r
#include "distributed_algebra_includes/MatrixGrandSum.f90"
#undef LMAT
#undef TEMP
    END IF
#undef MPIDATATYPE
  END SUBROUTINE MatrixGrandSum_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Elementwise multiplication. C_ij = A_ij * B_ij.
  !> Also known as a Hadamard product.
  RECURSIVE SUBROUTINE PairwiseMultiplyMatrix_ps(matA, matB, matC)
    !! Parameters
    !> Matrix A.
    TYPE(Matrix_ps), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_ps), INTENT(IN)  :: matB
    !> matC = MatA mult MatB.
    TYPE(Matrix_ps), INTENT(INOUT)  :: matC
    !! Local Data
    TYPE(Matrix_ps) :: converted_matrix
    INTEGER :: II, JJ

    IF (matA%is_complex .AND. .NOT. matB%is_complex) THEN
       CALL ConvertMatrixToComplex(matB, converted_matrix)
       CALL PairwiseMultiplyMatrix(matA, converted_matrix, matC)
       CALL DestructMatrix(converted_matrix)
    ELSE IF (.NOT. matA%is_complex .AND. matB%is_complex) THEN
       CALL ConvertMatrixToComplex(matA, converted_matrix)
       CALL PairwiseMultiplyMatrix(converted_matrix, matB, matC)
       CALL DestructMatrix(converted_matrix)
    ELSE IF (matA%is_complex .AND. matB%is_complex) THEN
#define LMAT local_data_c
#include "distributed_algebra_includes/PairwiseMultiply.f90"
#undef LMAT
    ELSE
#define LMAT local_data_r
#include "distributed_algebra_includes/PairwiseMultiply.f90"
#undef LMAT
    END IF
  END SUBROUTINE PairwiseMultiplyMatrix_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the norm of a distributed sparse matrix along the rows.
  FUNCTION MatrixNorm_ps(this) RESULT(norm_value)
    !! Parameters
    !> The matrix to compute the norm of.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The norm value of the full distributed sparse matrix.
    REAL(NTREAL) :: norm_value
    !! Local Data
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: local_norm
    TYPE(Matrix_lsr) :: merged_local_data_r
    TYPE(Matrix_lsc) :: merged_local_data_c
    INTEGER :: ierr

    IF (this%is_complex) THEN
#define LMAT merged_local_data_c
#include "distributed_algebra_includes/MatrixNorm.f90"
#undef LMAT
    ELSE
#define LMAT merged_local_data_r
#include "distributed_algebra_includes/MatrixNorm.f90"
#undef LMAT
    END IF
  END FUNCTION MatrixNorm_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> product = dot(Matrix A,Matrix B)
  !> Note that a dot product is the sum of elementwise multiplication, not
  !> traditional matrix multiplication.
  SUBROUTINE DotMatrix_psr(matA, matB, product)
    !> Matrix A.
    TYPE(Matrix_ps), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_ps), INTENT(IN)  :: matB
    !> The dot product.
    REAL(NTREAL), INTENT(OUT) :: product

    INCLUDE "distributed_algebra_includes/DotMatrix.f90"
  END SUBROUTINE DotMatrix_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> product = dot(Matrix A,Matrix B)
  !> Note that a dot product is the sum of elementwise multiplication, not
  !> traditional matrix multiplication.
  SUBROUTINE DotMatrix_psc(matA, matB, product)
    !> Matrix A.
    TYPE(Matrix_ps), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_ps), INTENT(IN)  :: matB
    !> The dot product.
    COMPLEX(NTCOMPLEX), INTENT(OUT) :: product

    INCLUDE "distributed_algebra_includes/DotMatrix.f90"
  END SUBROUTINE DotMatrix_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Matrix B = alpha*Matrix A + Matrix B (AXPY)
  !> This will utilize the sparse vector increment routine.
  RECURSIVE SUBROUTINE IncrementMatrix_ps(matA, matB, alpha_in, threshold_in)
    !> Matrix A.
    TYPE(Matrix_ps), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_ps), INTENT(INOUT)  :: matB
    !> Multiplier (default= 1.0).
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    !> For flushing values to zero (default=0).
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    !! Local Data
    TYPE(Matrix_ps) :: converted_matrix
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: threshold
    INTEGER :: II, JJ

    !! Optional Parameters
    IF (.NOT. PRESENT(alpha_in)) THEN
       alpha = 1.0d+0
    ELSE
       alpha = alpha_in
    END IF
    IF (.NOT. PRESENT(threshold_in)) THEN
       threshold = 0
    ELSE
       threshold = threshold_in
    END IF

    IF (matA%is_complex .AND. .NOT. matB%is_complex) THEN
       CALL ConvertMatrixToComplex(matB, converted_matrix)
       CALL IncrementMatrix(matA, converted_matrix, alpha, threshold)
       CALL CopyMatrix(converted_matrix, matB)
    ELSE IF (.NOT. matA%is_complex .AND. matB%is_complex) THEN
       CALL ConvertMatrixToComplex(matA, converted_matrix)
       CALL IncrementMatrix(converted_matrix, matB, alpha, threshold)
    ELSE IF (matA%is_complex .AND. matB%is_complex) THEN
#define LMAT local_data_c
#include "distributed_algebra_includes/IncrementMatrix.f90"
#undef LMAT
    ELSE
#define LMAT local_data_r
#include "distributed_algebra_includes/IncrementMatrix.f90"
#undef LMAT
    END IF

    CALL DestructMatrix(converted_matrix)

  END SUBROUTINE IncrementMatrix_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a distributed sparse matrix by a constant.
  SUBROUTINE ScaleMatrix_psr(this, constant)
    !> Matrix to scale.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> A constant scale factor.
    REAL(NTREAL), INTENT(IN) :: constant
    !! Local Data
    INTEGER :: II, JJ

    IF (this%is_complex) THEN
#define LOCALDATA local_data_c
#include "distributed_algebra_includes/ScaleMatrix.f90"
#undef LOCALDATA
    ELSE
#define LOCALDATA local_data_r
#include "distributed_algebra_includes/ScaleMatrix.f90"
#undef LOCALDATA
    END IF

  END SUBROUTINE ScaleMatrix_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a distributed sparse matrix by a constant.
  RECURSIVE SUBROUTINE ScaleMatrix_psc(this, constant)
    !> Matrix to scale.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> A constant scale factor.
    COMPLEX(NTCOMPLEX), INTENT(IN) :: constant
    !! Local Data
    TYPE(Matrix_ps) :: this_c
    INTEGER :: II, JJ

    IF (this%is_complex) THEN
#define LOCALDATA local_data_c
#include "distributed_algebra_includes/ScaleMatrix.f90"
#undef LOCALDATA
    ELSE
       CALL ConvertMatrixToComplex(this, this_c)
       CALL ScaleMatrix_psc(this_c, constant)
       CALL CopyMatrix(this_c, this)
       CALL DestructMatrix(this_c)
    END IF

  END SUBROUTINE ScaleMatrix_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the trace of the matrix.
  SUBROUTINE MatrixTrace_psr(this, trace_value)
    !! Parameters
    !> The matrix to compute the trace of.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The trace value of the full distributed sparse matrix.
    REAL(NTREAL), INTENT(OUT) :: trace_value
    !! Local data
    TYPE(TripletList_r) :: triplet_list_r
    TYPE(TripletList_c) :: triplet_list_c
    !! Counters/Temporary
    INTEGER :: counter
    TYPE(Matrix_lsr) :: merged_local_data_r
    TYPE(Matrix_lsc) :: merged_local_data_c
    INTEGER :: ierr

    IF (this%is_complex) THEN
#define TLIST triplet_list_c
#define LMAT merged_local_data_c
#define MPIDATATYPE MPINTCOMPLEX
#include "distributed_algebra_includes/MatrixTrace.f90"
#undef MPIDATATYPE
#undef LMAT
#undef TLIST
    ELSE
#define TLIST triplet_list_r
#define LMAT merged_local_data_r
#define MPIDATATYPE MPINTREAL
#include "distributed_algebra_includes/MatrixTrace.f90"
#undef MPIDATATYPE
#undef LMAT
#undef TLIST
    END IF
  END SUBROUTINE MatrixTrace_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE PSMatrixAlgebraModule
