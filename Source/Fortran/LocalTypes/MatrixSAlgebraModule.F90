! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !> A module for performing linear algebra using sparse matrices.
! MODULE MatrixSAlgebraModule
!   USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
!   USE MatrixDModule, ONLY : Matrix_ldr, Matrix_ldc, ConstructMatrixDFromS, &
!        & ConstructMatrixSFromD, MultiplyMatrix, DestructMatrix
!   USE MatrixMemoryPoolModule, ONLY : MatrixMemoryPool_lr, MatrixMemoryPool_lc, &
!        & DestructMatrixMemoryPool, CheckMemoryPoolValidity, SetPoolSparsity, &
!        & ConstructMatrixMemoryPool
!   USE MatrixSModule, ONLY: Matrix_lsr, Matrix_lsc, DestructMatrix, CopyMatrix, &
!        & TransposeMatrix, PrintMatrix, ConstructMatrixFromTripletList, &
!        & ConstructEmptyMatrix
!   USE VectorSModule, ONLY : AddSparseVectors, DotSparseVectors, &
!        & PairwiseMultiplyVectors
!   USE TripletListModule, ONLY: TripletList_r, TripletList_c, SortTripletList, &
!        & DestructTripletList, ConstructTripletList
!   USE TimerModule, ONLY : StartTimer, StopTimer
!   IMPLICIT NONE
!   PRIVATE
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   PUBLIC :: ScaleMatrix
!   PUBLIC :: IncrementMatrix
!   PUBLIC :: DotMatrix
!   PUBLIC :: PairwiseMultiplyMatrix
!   PUBLIC :: MatrixMultiply
!   PUBLIC :: MatrixColumnNorm
!   PUBLIC :: MatrixNorm
!   PUBLIC :: MatrixGrandSum
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   INTERFACE ScaleMatrix
!      MODULE PROCEDURE ScaleMatrix_lsr
!      MODULE PROCEDURE ScaleMatrix_lsc
!   END INTERFACE
!   INTERFACE IncrementMatrix
!      MODULE PROCEDURE IncrementMatrix_lsr
!      MODULE PROCEDURE IncrementMatrix_lsc
!   END INTERFACE
!   INTERFACE DotMatrix
!      MODULE PROCEDURE DotMatrix_lsr
!      MODULE PROCEDURE DotMatrix_lsc
!   END INTERFACE
!   INTERFACE PairwiseMultiplyMatrix
!      MODULE PROCEDURE PairwiseMultiplyMatrix_lsr
!      MODULE PROCEDURE PairwiseMultiplyMatrix_lsc
!   END INTERFACE
!   INTERFACE MatrixMultiply
!      MODULE PROCEDURE GemmMatrix_lsr
!      MODULE PROCEDURE GemmMatrix_lsc
!   END INTERFACE
!   INTERFACE MatrixColumnNorm
!      MODULE PROCEDURE MatrixColumnNorm_lsr
!      MODULE PROCEDURE MatrixColumnNorm_lsc
!   END INTERFACE
!   INTERFACE MatrixNorm
!      MODULE PROCEDURE MatrixNorm_lsr
!      MODULE PROCEDURE MatrixNorm_lsc
!   END INTERFACE
!   INTERFACE MatrixGrandSum
!      MODULE PROCEDURE MatrixGrandSum_lsr
!      MODULE PROCEDURE MatrixGrandSum_lsc
!   END INTERFACE
!   INTERFACE MultiplyBlock
!      MODULE PROCEDURE MultiplyBlock_lsr
!      MODULE PROCEDURE MultiplyBlock_lsc
!   END INTERFACE
!   INTERFACE PruneList
!      MODULE PROCEDURE PruneList_lsr
!      MODULE PROCEDURE PruneList_lsc
!   END INTERFACE
!   INTERFACE SparseBranch
!      MODULE PROCEDURE SparseBranch_lsr
!      MODULE PROCEDURE SparseBranch_lsc
!   END INTERFACE
!   INTERFACE DenseBranch
!      MODULE PROCEDURE DenseBranch_lsr
!      MODULE PROCEDURE DenseBranch_lsc
!   END INTERFACE
! CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Will scale a sparse matrix by a constant.
!   !! @param[inout] matA Matrix A.
!   !! @param[in] constant scale factor.
!   PURE SUBROUTINE ScaleMatrix_lsr(matA,constant)
!     !! Parameters
!     TYPE(Matrix_lsr), INTENT(INOUT) :: matA
!     REAL(NTREAL), INTENT(IN) :: constant
!
!     INCLUDE "sparse_includes/ScaleMatrix.f90"
!   END SUBROUTINE ScaleMatrix_lsr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Matrix B = alpha*Matrix A + Matrix B (AXPY).
!   !! This will utilize the sparse vector addition routine.
!   !! @param[in] matA Matrix A.
!   !! @param[in,out] matB Matrix B.
!   !! @param[in] alpha_in multiplier (optional, default=1.0)
!   !! @param[in] threshold_in for flushing values to zero. (Optional, default=0).
!   PURE SUBROUTINE IncrementMatrix_lsr(matA, matB, alpha_in, threshold_in)
!     !! Parameters
!     TYPE(Matrix_lsr), INTENT(IN)  :: matA
!     TYPE(Matrix_lsr), INTENT(INOUT) :: matB
!     REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
!     REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
!     !! Local Variables
!     TYPE(Matrix_lsr) :: matC
!
!     INCLUDE "sparse_includes/IncrementMatrix.f90"
!   END SUBROUTINE IncrementMatrix_lsr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Pairwise Multiply two matrices.
!   !! This will utilize the sparse vector pairwise routine.
!   !! @param[in] matA Matrix A.
!   !! @param[in] matB Matrix B.
!   !! @param[in,out] matC = MatA mult MatB.
!   PURE SUBROUTINE PairwiseMultiplyMatrix_lsr(matA, matB, matC)
!     !! Parameters
!     TYPE(Matrix_lsr), INTENT(IN)  :: matA
!     TYPE(Matrix_lsr), INTENT(IN) :: matB
!     TYPE(Matrix_lsr), INTENT(INOUT) :: matC
!     !! Local Variables
!     TYPE(Matrix_lsr) :: TempMat
!
!     INCLUDE "sparse_includes/PairwiseMultiplyMatrix.f90"
!   END SUBROUTINE PairwiseMultiplyMatrix_lsr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Product = sum(MatA[ij]*MatB[ij])
!   !! @param[in] matA Matrix A.
!   !! @param[in] matB Matrix B.
!   !! @result product
!   PURE FUNCTION DotMatrix_lsr(matA, matB) RESULT(product)
!     !! Parameters
!     TYPE(Matrix_lsr), INTENT(IN) :: matA
!     TYPE(Matrix_lsr), INTENT(IN) :: matB
!     REAL(NTREAL) :: product
!     !! Local Variables
!     TYPE(Matrix_lsr) :: matC
!
!     INCLUDE "sparse_includes/DotMatrix.f90"
!   END FUNCTION DotMatrix_lsr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Multiply two matrices together, and add to the third.
!   !! C := alpha*matA*op( matB ) + beta*matC
!   !! @param[in] matA Matrix A.
!   !! @param[in] matB Matrix B.
!   !! @param[out] matC = alpha*matA*op( matB ) + beta*matC.
!   !! @param[in] IsATransposed_in true if A is already transposed.
!   !! @param[in] IsBTransposed_in true if B is already transposed.
!   !! @param[in] alpha_in scales the multiplication.
!   !! @param[in] beta_in scales matrix we sum on to.
!   !! @param[in] threshold_in for flushing values to zero. Default value is 0.0.
!   !! @param[inout] blocked_memory_pool_in an optional memory pool for doing the
!   !! calculation.
!   SUBROUTINE GemmMatrix_lsr(matA, matB, matC, IsATransposed_in, IsBTransposed_in, &
!        & alpha_in, beta_in, threshold_in, blocked_memory_pool_in)
!     !! Parameters
!     TYPE(Matrix_lsr), INTENT(IN)  :: matA
!     TYPE(Matrix_lsr), INTENT(IN)  :: matB
!     TYPE(Matrix_lsr), INTENT(INOUT) :: matC
!     LOGICAL, OPTIONAL, INTENT(IN) :: IsATransposed_in
!     LOGICAL, OPTIONAL, INTENT(IN) :: IsBTransposed_in
!     REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
!     REAL(NTREAL), OPTIONAL, INTENT(IN) :: beta_in
!     REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
!     TYPE(MatrixMemoryPool_lr), OPTIONAL, &
!          & INTENT(INOUT), TARGET :: blocked_memory_pool_in
!     !! Intermediate Data
!     TYPE(Matrix_lsr) :: matAB
!     LOGICAL :: IsATransposed, IsBTransposed
!     REAL(NTREAL) :: alpha
!     REAL(NTREAL) :: beta
!     REAL(NTREAL) :: threshold
!     TYPE(MatrixMemoryPool_lr) :: blocked_memory_pool
!
!     INCLUDE "sparse_includes/GemmMatrix.f90"
!   END SUBROUTINE GemmMatrix_lsr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Compute the norm of a sparse matrix along the columns.
!   !! @param[in] this the matrix to compute the norm of.
!   !! @param[out] norm_per_column the norm value for each column in this matrix.
!   PURE SUBROUTINE MatrixColumnNorm_lsr(this, norm_per_column)
!     !! Parameters
!     TYPE(Matrix_lsr), INTENT(IN) :: this
!     REAL(NTREAL), DIMENSION(this%columns), INTENT(OUT) :: norm_per_column
!     !! Local Data
!     REAL(NTREAL) :: temp_value
!
!     INCLUDE "sparse_includes/MatrixColumnNorm.f90"
!   END SUBROUTINE MatrixColumnNorm_lsr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Compute the 1 norm of a sparse matrix.
!   !! @param[in] this the matrix to compute the norm of.
!   !! @result norm the matrix.
!   PURE FUNCTION MatrixNorm_lsr(this) RESULT(norm)
!     !! Parameters
!     TYPE(Matrix_lsr), INTENT(IN) :: this
!     REAL(NTREAL) :: norm
!     !! Local Variables
!     REAL(NTREAL), DIMENSION(this%columns) :: column
!
!     INCLUDE "sparse_includes/MatrixNorm.f90"
!
!   END FUNCTION MatrixNorm_lsr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Sum the elements of a matrix
!   !! @param[in] this the matrix to sum
!   !! @result sum_value the sum of the matrix elements
!   PURE FUNCTION MatrixGrandSum_lsr(this) RESULT(sum_value)
!     !! Parameters
!     TYPE(Matrix_lsr), INTENT(IN) :: this
!     REAL(NTREAL) :: sum_value
!
!     INCLUDE "sparse_includes/MatrixGrandSum.f90"
!
!   END FUNCTION MatrixGrandSum_lsr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   PURE SUBROUTINE SparseBranch_lsr(matA, matB, matC, IsATransposed, IsBTransposed, &
!        & alpha, threshold, blocked_memory_pool)
!     !! Parameters
!     TYPE(Matrix_lsr), INTENT(IN)  :: matA
!     TYPE(Matrix_lsr), INTENT(IN)  :: matB
!     TYPE(Matrix_lsr), INTENT(INOUT) :: matC
!     LOGICAL, INTENT(IN) :: IsATransposed
!     LOGICAL, INTENT(IN) :: IsBTransposed
!     REAL(NTREAL), INTENT(IN) :: alpha
!     REAL(NTREAL), INTENT(IN) :: threshold
!     TYPE(MatrixMemoryPool_lr), INTENT(INOUT) :: blocked_memory_pool
!     !! Local Data
!     TYPE(Matrix_lsr) :: matAT, matBT
!
!     INCLUDE "sparse_includes/SparseBranch.f90"
!   END SUBROUTINE SparseBranch_lsr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   SUBROUTINE DenseBranch_lsr(matA, matB, matC, IsATransposed, IsBTransposed, &
!        & alpha, threshold)
!     !! Parameters
!     TYPE(Matrix_lsr), INTENT(IN)  :: matA
!     TYPE(Matrix_lsr), INTENT(IN)  :: matB
!     TYPE(Matrix_lsr), INTENT(INOUT) :: matC
!     LOGICAL, INTENT(IN) :: IsATransposed
!     LOGICAL, INTENT(IN) :: IsBTransposed
!     REAL(NTREAL), INTENT(IN) :: alpha
!     REAL(NTREAL), INTENT(IN) :: threshold
!     !! Local Data
!     TYPE(Matrix_lsr) :: untransposedMatA
!     TYPE(Matrix_lsr) :: untransposedMatB
!     TYPE(Matrix_ldr) :: DenseA
!     TYPE(Matrix_ldr) :: DenseB
!     TYPE(Matrix_ldr) :: DenseC
!
!     INCLUDE "sparse_includes/DenseBranch.f90"
!   END SUBROUTINE DenseBranch_lsr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   PURE SUBROUTINE MultiplyBlock_lsr(matAT,matBT,memorypool)
!     !! Parameters
!     TYPE(Matrix_lsr), INTENT(IN)  :: matAT
!     TYPE(Matrix_lsr), INTENT(IN)  :: matBT
!     TYPE(MatrixMemoryPool_lr), INTENT(INOUT) :: memorypool
!     !! Temp Variables
!     REAL(NTREAL) :: temp_value_a, temp_value_b, temp_value_c
!
!     INCLUDE "sparse_includes/MultiplyBlock.f90"
!   END SUBROUTINE MultiplyBlock_lsr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   PURE SUBROUTINE PruneList_lsr(memorypool,alpha,threshold, &
!        & mat_c_columns, mat_c_rows, matAB)
!     !! Parameters
!     TYPE(MatrixMemoryPool_lr), INTENT(INOUT) :: memorypool
!     REAL(NTREAL), INTENT(IN) :: alpha
!     REAL(NTREAL), INTENT(IN) :: threshold
!     INTEGER, INTENT(IN) :: mat_c_columns
!     INTEGER, INTENT(IN) :: mat_c_rows
!     TYPE(Matrix_lsr), INTENT(INOUT) :: matAB
!     !! Local data
!     REAL(NTREAL) :: working_value
!     TYPE(TripletList_r) :: unsorted_pruned_list
!     TYPE(TripletList_r) :: sorted_pruned_list
!
!     INCLUDE "sparse_includes/PruneList.f90"
!   END SUBROUTINE PruneList_lsr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Will scale a sparse matrix by a constant.
!   !! @param[inout] matA Matrix A.
!   !! @param[in] constant scale factor.
!   PURE SUBROUTINE ScaleMatrix_lsc(matA,constant)
!     !! Parameters
!     TYPE(Matrix_lsc), INTENT(INOUT) :: matA
!     REAL(NTREAL), INTENT(IN) :: constant
!
!     INCLUDE "sparse_includes/ScaleMatrix.f90"
!   END SUBROUTINE ScaleMatrix_lsc
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Matrix B = alpha*Matrix A + Matrix B (AXPY).
!   !! This will utilize the sparse vector addition routine.
!   !! @param[in] matA Matrix A.
!   !! @param[in,out] matB Matrix B.
!   !! @param[in] alpha_in multiplier (optional, default=1.0)
!   !! @param[in] threshold_in for flushing values to zero. (Optional, default=0).
!   PURE SUBROUTINE IncrementMatrix_lsc(matA, matB, alpha_in, threshold_in)
!     !! Parameters
!     TYPE(Matrix_lsc), INTENT(IN)  :: matA
!     TYPE(Matrix_lsc), INTENT(INOUT) :: matB
!     REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
!     REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
!     !! Local Variables
!     TYPE(Matrix_lsc) :: matC
!
!     INCLUDE "sparse_includes/IncrementMatrix.f90"
!   END SUBROUTINE IncrementMatrix_lsc
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Pairwise Multiply two matrices.
!   !! This will utilize the sparse vector pairwise routine.
!   !! @param[in] matA Matrix A.
!   !! @param[in] matB Matrix B.
!   !! @param[in,out] matC = MatA mult MatB.
!   PURE SUBROUTINE PairwiseMultiplyMatrix_lsc(matA, matB, matC)
!     !! Parameters
!     TYPE(Matrix_lsc), INTENT(IN)  :: matA
!     TYPE(Matrix_lsc), INTENT(IN) :: matB
!     TYPE(Matrix_lsc), INTENT(INOUT) :: matC
!     !! Local Variables
!     TYPE(Matrix_lsc) :: TempMat
!
!     INCLUDE "sparse_includes/PairwiseMultiplyMatrix.f90"
!   END SUBROUTINE PairwiseMultiplyMatrix_lsc
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Product = sum(MatA[ij]*MatB[ij])
!   !! @param[in] matA Matrix A.
!   !! @param[in] matB Matrix B.
!   !! @result product
!   PURE FUNCTION DotMatrix_lsc(matA, matB) RESULT(product)
!     !! Parameters
!     TYPE(Matrix_lsc), INTENT(IN) :: matA
!     TYPE(Matrix_lsc), INTENT(IN) :: matB
!     COMPLEX(NTCOMPLEX) :: product
!     !! Local Variables
!     TYPE(Matrix_lsc) :: matC
!
!     INCLUDE "sparse_includes/DotMatrix.f90"
!   END FUNCTION DotMatrix_lsc
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Multiply two matrices together, and add to the third.
!   !! C := alpha*matA*op( matB ) + beta*matC
!   !! @param[in] matA Matrix A.
!   !! @param[in] matB Matrix B.
!   !! @param[out] matC = alpha*matA*op( matB ) + beta*matC.
!   !! @param[in] IsATransposed_in true if A is already transposed.
!   !! @param[in] IsBTransposed_in true if B is already transposed.
!   !! @param[in] alpha_in scales the multiplication.
!   !! @param[in] beta_in scales matrix we sum on to.
!   !! @param[in] threshold_in for flushing values to zero. Default value is 0.0.
!   !! @param[inout] blocked_memory_pool_in an optional memory pool for doing the
!   !! calculation.
!   SUBROUTINE GemmMatrix_lsc(matA, matB, matC, IsATransposed_in, IsBTransposed_in, &
!        & alpha_in, beta_in, threshold_in, blocked_memory_pool_in)
!     !! Parameters
!     TYPE(Matrix_lsc), INTENT(IN)  :: matA
!     TYPE(Matrix_lsc), INTENT(IN)  :: matB
!     TYPE(Matrix_lsc), INTENT(INOUT) :: matC
!     LOGICAL, OPTIONAL, INTENT(IN) :: IsATransposed_in
!     LOGICAL, OPTIONAL, INTENT(IN) :: IsBTransposed_in
!     REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
!     REAL(NTREAL), OPTIONAL, INTENT(IN) :: beta_in
!     REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
!     TYPE(MatrixMemoryPool_lc), OPTIONAL, &
!          & INTENT(INOUT), TARGET :: blocked_memory_pool_in
!     !! Intermediate Data
!     TYPE(Matrix_lsc) :: matAB
!     LOGICAL :: IsATransposed, IsBTransposed
!     REAL(NTREAL) :: alpha
!     REAL(NTREAL) :: beta
!     REAL(NTREAL) :: threshold
!     TYPE(MatrixMemoryPool_lc) :: blocked_memory_pool
!
!     INCLUDE "sparse_includes/GemmMatrix.f90"
!   END SUBROUTINE GemmMatrix_lsc
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Compute the norm of a sparse matrix along the columns.
!   !! @param[in] this the matrix to compute the norm of.
!   !! @param[out] norm_per_column the norm value for each column in this matrix.
!   PURE SUBROUTINE MatrixColumnNorm_lsc(this, norm_per_column)
!     !! Parameters
!     TYPE(Matrix_lsc), INTENT(IN) :: this
!     REAL(NTREAL), DIMENSION(this%columns), INTENT(OUT) :: norm_per_column
!     !! Local Data
!     COMPLEX(NTCOMPLEX)  :: temp_value
!
!     INCLUDE "sparse_includes/MatrixColumnNorm.f90"
!   END SUBROUTINE MatrixColumnNorm_lsc
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Compute the 1 norm of a sparse matrix.
!   !! @param[in] this the matrix to compute the norm of.
!   !! @result norm the matrix.
!   PURE FUNCTION MatrixNorm_lsc(this) RESULT(norm)
!     !! Parameters
!     TYPE(Matrix_lsc), INTENT(IN) :: this
!     REAL(NTREAL) :: norm
!     !! Local Variables
!     REAL(NTREAL), DIMENSION(this%columns) :: column
!
!     INCLUDE "sparse_includes/MatrixNorm.f90"
!
!   END FUNCTION MatrixNorm_lsc
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Sum the elements of a matrix
!   !! @param[in] this the matrix to sum
!   !! @result sum_value the sum of the matrix elements
!   PURE FUNCTION MatrixGrandSum_lsc(this) RESULT(sum_value)
!     !! Parameters
!     TYPE(Matrix_lsc), INTENT(IN) :: this
!     COMPLEX(NTCOMPLEX) :: sum_value
!
!     INCLUDE "sparse_includes/MatrixGrandSum.f90"
!
!   END FUNCTION MatrixGrandSum_lsc
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   PURE SUBROUTINE SparseBranch_lsc(matA, matB, matC, IsATransposed, IsBTransposed, &
!        & alpha, threshold, blocked_memory_pool)
!     !! Parameters
!     TYPE(Matrix_lsc), INTENT(IN)  :: matA
!     TYPE(Matrix_lsc), INTENT(IN)  :: matB
!     TYPE(Matrix_lsc), INTENT(INOUT) :: matC
!     LOGICAL, INTENT(IN) :: IsATransposed
!     LOGICAL, INTENT(IN) :: IsBTransposed
!     REAL(NTREAL), INTENT(IN) :: alpha
!     REAL(NTREAL), INTENT(IN) :: threshold
!     TYPE(MatrixMemoryPool_lc), INTENT(INOUT) :: blocked_memory_pool
!     !! Local Data
!     TYPE(Matrix_lsc) :: matAT, matBT
!
!     INCLUDE "sparse_includes/SparseBranch.f90"
!   END SUBROUTINE SparseBranch_lsc
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   SUBROUTINE DenseBranch_lsc(matA, matB, matC, IsATransposed, IsBTransposed, &
!        & alpha, threshold)
!     !! Parameters
!     TYPE(Matrix_lsc), INTENT(IN)  :: matA
!     TYPE(Matrix_lsc), INTENT(IN)  :: matB
!     TYPE(Matrix_lsc), INTENT(INOUT) :: matC
!     LOGICAL, INTENT(IN) :: IsATransposed
!     LOGICAL, INTENT(IN) :: IsBTransposed
!     REAL(NTREAL), INTENT(IN) :: alpha
!     REAL(NTREAL), INTENT(IN) :: threshold
!     !! Local Data
!     TYPE(Matrix_lsc) :: untransposedMatA
!     TYPE(Matrix_lsc) :: untransposedMatB
!     TYPE(Matrix_ldc) :: DenseA
!     TYPE(Matrix_ldc) :: DenseB
!     TYPE(Matrix_ldc) :: DenseC
!
!     INCLUDE "sparse_includes/DenseBranch.f90"
!   END SUBROUTINE DenseBranch_lsc
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   PURE SUBROUTINE MultiplyBlock_lsc(matAT,matBT,memorypool)
!     !! Parameters
!     TYPE(Matrix_lsc), INTENT(IN)  :: matAT
!     TYPE(Matrix_lsc), INTENT(IN)  :: matBT
!     TYPE(MatrixMemoryPool_lc), INTENT(INOUT) :: memorypool
!     !! Temp Variables
!     COMPLEX(NTCOMPLEX) :: temp_value_a, temp_value_b, temp_value_c
!
!     INCLUDE "sparse_includes/MultiplyBlock.f90"
!   END SUBROUTINE MultiplyBlock_lsc
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   PURE SUBROUTINE PruneList_lsc(memorypool,alpha,threshold, &
!        & mat_c_columns, mat_c_rows, matAB)
!     !! Parameters
!     TYPE(MatrixMemoryPool_lc), INTENT(INOUT) :: memorypool
!     REAL(NTREAL), INTENT(IN) :: alpha
!     REAL(NTREAL), INTENT(IN) :: threshold
!     INTEGER, INTENT(IN) :: mat_c_columns
!     INTEGER, INTENT(IN) :: mat_c_rows
!     TYPE(Matrix_lsc), INTENT(INOUT) :: matAB
!     !! Local data
!     COMPLEX(NTCOMPLEX) :: working_value
!     TYPE(TripletList_c) :: unsorted_pruned_list
!     TYPE(TripletList_c) :: sorted_pruned_list
!
!     INCLUDE "sparse_includes/PruneList.f90"
!   END SUBROUTINE PruneList_lsc
! END MODULE MatrixSAlgebraModule
