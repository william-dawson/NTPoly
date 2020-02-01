!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for performing linear algebra using sparse matrices.
MODULE SMatrixAlgebraModule
  USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
  USE DMatrixModule, ONLY : Matrix_ldr, Matrix_ldc, ConstructMatrixDFromS, &
       & ConstructMatrixSFromD, CopyMatrix, MultiplyMatrix, TransposeMatrix, &
       & DestructMatrix
  USE MatrixMemoryPoolModule, ONLY : MatrixMemoryPool_lr, MatrixMemoryPool_lc, &
       & DestructMatrixMemoryPool, CheckMemoryPoolValidity, SetPoolSparsity, &
       & ConstructMatrixMemoryPool
  USE SMatrixModule, ONLY: Matrix_lsr, Matrix_lsc, DestructMatrix, CopyMatrix, &
       & TransposeMatrix, ConjugateMatrix, ConstructMatrixFromTripletList, &
       & ConstructEmptyMatrix
  USE SVectorModule, ONLY : AddSparseVectors, PairwiseMultiplyVectors
  USE TripletListModule, ONLY: TripletList_r, TripletList_c, SortTripletList, &
       & DestructTripletList, ConstructTripletList
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ScaleMatrix
  PUBLIC :: IncrementMatrix
  PUBLIC :: DotMatrix
  PUBLIC :: PairwiseMultiplyMatrix
  PUBLIC :: MatrixMultiply
  PUBLIC :: MatrixColumnNorm
  PUBLIC :: MatrixNorm
  PUBLIC :: MatrixGrandSum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE ScaleMatrix
     MODULE PROCEDURE ScaleMatrix_lsr
     MODULE PROCEDURE ScaleMatrix_lsc
     MODULE PROCEDURE ScaleMatrix_lsc_c
  END INTERFACE
  INTERFACE IncrementMatrix
     MODULE PROCEDURE IncrementMatrix_lsr
     MODULE PROCEDURE IncrementMatrix_lsc
  END INTERFACE
  INTERFACE DotMatrix
     MODULE PROCEDURE DotMatrix_lsr
     MODULE PROCEDURE DotMatrix_lsc
  END INTERFACE
  INTERFACE PairwiseMultiplyMatrix
     MODULE PROCEDURE PairwiseMultiplyMatrix_lsr
     MODULE PROCEDURE PairwiseMultiplyMatrix_lsc
  END INTERFACE
  INTERFACE MatrixMultiply
     MODULE PROCEDURE GemmMatrix_lsr
     MODULE PROCEDURE GemmMatrix_lsc
  END INTERFACE
  INTERFACE MatrixColumnNorm
     MODULE PROCEDURE MatrixColumnNorm_lsr
     MODULE PROCEDURE MatrixColumnNorm_lsc
  END INTERFACE
  INTERFACE MatrixNorm
     MODULE PROCEDURE MatrixNorm_lsr
     MODULE PROCEDURE MatrixNorm_lsc
  END INTERFACE
  INTERFACE MatrixGrandSum
     MODULE PROCEDURE MatrixGrandSum_lsr
     MODULE PROCEDURE MatrixGrandSum_lsc
  END INTERFACE
  INTERFACE MultiplyBlock
     MODULE PROCEDURE MultiplyBlock_lsr
     MODULE PROCEDURE MultiplyBlock_lsc
  END INTERFACE
  INTERFACE PruneList
     MODULE PROCEDURE PruneList_lsr
     MODULE PROCEDURE PruneList_lsc
  END INTERFACE
  INTERFACE SparseBranch
     MODULE PROCEDURE SparseBranch_lsr
     MODULE PROCEDURE SparseBranch_lsc
  END INTERFACE
  INTERFACE DenseBranch
     MODULE PROCEDURE DenseBranch_lsr
     MODULE PROCEDURE DenseBranch_lsc
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a sparse matrix by a constant.
  PURE SUBROUTINE ScaleMatrix_lsr(matA,constant)
    !> Matrix A.
    TYPE(Matrix_lsr), INTENT(INOUT) :: matA
    !> Constant scale factor.
    REAL(NTREAL), INTENT(IN) :: constant

    INCLUDE "sparse_includes/ScaleMatrix.f90"
  END SUBROUTINE ScaleMatrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a sparse matrix by a constant.
  PURE SUBROUTINE ScaleMatrix_lsc(matA,constant)
    !> Matrix A.
    TYPE(Matrix_lsc), INTENT(INOUT) :: matA
    !> Constant scale factor.
    REAL(NTREAL), INTENT(IN) :: constant

    INCLUDE "sparse_includes/ScaleMatrix.f90"
  END SUBROUTINE ScaleMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a sparse matrix by a constant.
  PURE SUBROUTINE ScaleMatrix_lsc_c(matA,constant)
    !> Matrix A.
    TYPE(Matrix_lsc), INTENT(INOUT) :: matA
    !> Constant scale factor.
    COMPLEX(NTCOMPLEX), INTENT(IN) :: constant

    INCLUDE "sparse_includes/ScaleMatrix.f90"
  END SUBROUTINE ScaleMatrix_lsc_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Matrix B = alpha*Matrix A + Matrix B (AXPY).
  !> This will utilize the sparse vector addition routine.
  PURE SUBROUTINE IncrementMatrix_lsr(matA, matB, alpha_in, threshold_in)
    !> Matrix A.
    TYPE(Matrix_lsr), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_lsr), INTENT(INOUT) :: matB
    !> Multiplier (default=1.0).
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    !> For flushing values to zero (default=0).
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    !! Local Variables
    TYPE(Matrix_lsr) :: matC

    INCLUDE "sparse_includes/IncrementMatrix.f90"
  END SUBROUTINE IncrementMatrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Matrix B = alpha*Matrix A + Matrix B (AXPY).
  !> This will utilize the sparse vector addition routine.
  PURE SUBROUTINE IncrementMatrix_lsc(matA, matB, alpha_in, threshold_in)
    !> Matrix A.
    TYPE(Matrix_lsc), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_lsc), INTENT(INOUT) :: matB
    !> Multiplier (default=1.0).
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    !> For flushing values to zero (default=0).
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    !! Local Variables
    TYPE(Matrix_lsc) :: matC

    INCLUDE "sparse_includes/IncrementMatrix.f90"
  END SUBROUTINE IncrementMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Pairwise Multiply two matrices.
  !> This will utilize the sparse vector pairwise multiply routine.
  PURE SUBROUTINE PairwiseMultiplyMatrix_lsr(matA, matB, matC)
    !> Matrix A.
    TYPE(Matrix_lsr), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_lsr), INTENT(IN) :: matB
    !> matC = MatA mult MatB.
    TYPE(Matrix_lsr), INTENT(INOUT) :: matC
    !! Local Variables
    TYPE(Matrix_lsr) :: TempMat

    INCLUDE "sparse_includes/PairwiseMultiplyMatrix.f90"
  END SUBROUTINE PairwiseMultiplyMatrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Pairwise Multiply two matrices.
  !> This will utilize the sparse vector pairwise routine.
  PURE SUBROUTINE PairwiseMultiplyMatrix_lsc(matA, matB, matC)
    !> Matrix A.
    TYPE(Matrix_lsc), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_lsc), INTENT(IN) :: matB
    !> matC = MatA mult MatB.
    TYPE(Matrix_lsc), INTENT(INOUT) :: matC
    !! Local Variables
    TYPE(Matrix_lsc) :: TempMat

    INCLUDE "sparse_includes/PairwiseMultiplyMatrix.f90"
  END SUBROUTINE PairwiseMultiplyMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Product = sum(MatA[ij]*MatB[ij])
  PURE SUBROUTINE DotMatrix_lsr(matA, matB, product)
    !> Matrix A.
    TYPE(Matrix_lsr), INTENT(IN) :: matA
    !> Matrix B.
    TYPE(Matrix_lsr), INTENT(IN) :: matB
    !> Dot product.
    REAL(NTREAL), INTENT(OUT) :: product
    !! Local Variables
    TYPE(Matrix_lsr) :: matC

    CALL PairwiseMultiplyMatrix(matA,matB,matC)

    CALL MatrixGrandSum(matC, product)
    CALL DestructMatrix(matC)

  END SUBROUTINE DotMatrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Product = sum(MatA^H[ij]*MatB[ij])
  PURE SUBROUTINE DotMatrix_lsc(matA, matB, product)
    !> Matrix A.
    TYPE(Matrix_lsc), INTENT(IN) :: matA
    !> Matrix B.
    TYPE(Matrix_lsc), INTENT(IN) :: matB
    !> Dot product.
    COMPLEX(NTCOMPLEX), INTENT(OUT) :: product
    !! Local Variables
    TYPE(Matrix_lsc) :: matC
    TYPE(Matrix_lsc) :: matAH

    CALL CopyMatrix(matA, matAH)
    CALL ConjugateMatrix(matAH)

    CALL PairwiseMultiplyMatrix(matAH, matB, matC)
    CALL MatrixGrandSum(matC, product)

    CALL DestructMatrix(matC)
    CALL DestructMatrix(matAH)

  END SUBROUTINE DotMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Multiply two matrices together, and add to the third.
  !> C := alpha*matA*op( matB ) + beta*matC
  SUBROUTINE GemmMatrix_lsr(matA, matB, matC, IsATransposed_in, &
       & IsBTransposed_in, alpha_in, beta_in, threshold_in, &
       & blocked_memory_pool_in)
    !> Matrix A.
    TYPE(Matrix_lsr), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_lsr), INTENT(IN)  :: matB
    !> matC = alpha*matA*op( matB ) + beta*matC.
    TYPE(Matrix_lsr), INTENT(INOUT) :: matC
    !> True if A is already transposed.
    LOGICAL, OPTIONAL, INTENT(IN) :: IsATransposed_in
    !> True if B is already transposed.
    LOGICAL, OPTIONAL, INTENT(IN) :: IsBTransposed_in
    !> Scales the multiplication.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    !> Scales matrix we sum on to.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: beta_in
    !> For flushing values to zero. Default value is 0.0.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    !> An optional memory pool for doing the calculation.
    TYPE(MatrixMemoryPool_lr), OPTIONAL, &
         & INTENT(INOUT), TARGET :: blocked_memory_pool_in
    !! Intermediate Data
    TYPE(Matrix_lsr) :: matAB
    LOGICAL :: IsATransposed, IsBTransposed
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: beta
    REAL(NTREAL) :: threshold
    TYPE(MatrixMemoryPool_lr) :: blocked_memory_pool

    INCLUDE "sparse_includes/GemmMatrix.f90"
  END SUBROUTINE GemmMatrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Multiply two matrices together, and add to the third.
  !> C := alpha*matA*op( matB ) + beta*matC
  SUBROUTINE GemmMatrix_lsc(matA, matB, matC, IsATransposed_in, &
       & IsBTransposed_in, alpha_in, beta_in, threshold_in, &
       & blocked_memory_pool_in)
    !> Matrix A.
    TYPE(Matrix_lsc), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_lsc), INTENT(IN)  :: matB
    !> matC = alpha*matA*op( matB ) + beta*matC.
    TYPE(Matrix_lsc), INTENT(INOUT) :: matC
    !> True if A is already transposed.
    LOGICAL, OPTIONAL, INTENT(IN) :: IsATransposed_in
    !> True if B is already transposed.
    LOGICAL, OPTIONAL, INTENT(IN) :: IsBTransposed_in
    !> Scales the multiplication.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    !> Scales matrix we sum on to.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: beta_in
    !> For flushing values to zero. Default value is 0.0.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    !> An optional memory pool for doing the calculation.
    TYPE(MatrixMemoryPool_lc), OPTIONAL, &
         & INTENT(INOUT), TARGET :: blocked_memory_pool_in
    !! Intermediate Data
    TYPE(Matrix_lsc) :: matAB
    LOGICAL :: IsATransposed, IsBTransposed
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: beta
    REAL(NTREAL) :: threshold
    TYPE(MatrixMemoryPool_lc) :: blocked_memory_pool

    INCLUDE "sparse_includes/GemmMatrix.f90"
  END SUBROUTINE GemmMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the norm of a sparse matrix along the columns.
  PURE SUBROUTINE MatrixColumnNorm_lsr(this, norm_per_column)
    !> The matrix to compute the norm of.
    TYPE(Matrix_lsr), INTENT(IN) :: this
    !> The norm value for each column in this matrix.
    REAL(NTREAL), DIMENSION(this%columns), INTENT(OUT) :: norm_per_column
    !! Local Data
    REAL(NTREAL) :: temp_value

    INCLUDE "sparse_includes/MatrixColumnNorm.f90"
  END SUBROUTINE MatrixColumnNorm_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the norm of a sparse matrix along the columns.
  PURE SUBROUTINE MatrixColumnNorm_lsc(this, norm_per_column)
    !> The matrix to compute the norm of.
    TYPE(Matrix_lsc), INTENT(IN) :: this
    !> The norm value for each column in this matrix.
    REAL(NTREAL), DIMENSION(this%columns), INTENT(OUT) :: norm_per_column
    !! Local Data
    COMPLEX(NTCOMPLEX)  :: temp_value

    INCLUDE "sparse_includes/MatrixColumnNorm.f90"
  END SUBROUTINE MatrixColumnNorm_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the 1 norm of a sparse matrix.
  PURE FUNCTION MatrixNorm_lsr(this) RESULT(norm)
    !> The matrix to compute the norm of.
    TYPE(Matrix_lsr), INTENT(IN) :: this
    !> The norm of the matrix.
    REAL(NTREAL) :: norm
    !! Local Variables
    REAL(NTREAL), DIMENSION(this%columns) :: column

    INCLUDE "sparse_includes/MatrixNorm.f90"

  END FUNCTION MatrixNorm_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the 1 norm of a sparse matrix.
  PURE FUNCTION MatrixNorm_lsc(this) RESULT(norm)
    !> The matrix to compute the norm of.
    TYPE(Matrix_lsc), INTENT(IN) :: this
    !> The norm of the matrix.
    REAL(NTREAL) :: norm
    !! Local Variables
    REAL(NTREAL), DIMENSION(this%columns) :: column

    INCLUDE "sparse_includes/MatrixNorm.f90"

  END FUNCTION MatrixNorm_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sum the elements of a matrix
  PURE SUBROUTINE MatrixGrandSum_lsr(this, sum_value)
    !> The matrix to sum
    TYPE(Matrix_lsr), INTENT(IN) :: this
    !> The sum of the matrix elements
    REAL(NTREAL), INTENT(OUT) :: sum_value

    INCLUDE "sparse_includes/MatrixGrandSum.f90"

  END SUBROUTINE MatrixGrandSum_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sum the elements of a matrix
  PURE SUBROUTINE MatrixGrandSum_lsc(this, sum_value)
    !> The matrix to sum
    TYPE(Matrix_lsc), INTENT(IN) :: this
    !> The sum of the matrix elements
    COMPLEX(NTCOMPLEX), INTENT(OUT) :: sum_value

    INCLUDE "sparse_includes/MatrixGrandSum.f90"

  END SUBROUTINE MatrixGrandSum_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Calculates the matrix product if using sparse-sparse algorithm.
  PURE SUBROUTINE SparseBranch_lsr(matA, matB, matC, IsATransposed, &
       & IsBTransposed, alpha, threshold, blocked_memory_pool)
    !> Matrix A.
    TYPE(Matrix_lsr), INTENT(IN)  :: matA
    !> Matrix B
    TYPE(Matrix_lsr), INTENT(IN)  :: matB
    !> matC = alpha*matA*op( matB ) + beta*matC.
    TYPE(Matrix_lsr), INTENT(INOUT) :: matC
    !> True if A is transposed.
    LOGICAL, INTENT(IN) :: IsATransposed
    !> True if B is transposed.
    LOGICAL, INTENT(IN) :: IsBTransposed
    !> Scaling value.
    REAL(NTREAL), INTENT(IN) :: alpha
    !> Threshold for flushing values.
    REAL(NTREAL), INTENT(IN) :: threshold
    !> Memory pool.
    TYPE(MatrixMemoryPool_lr), INTENT(INOUT) :: blocked_memory_pool
    !! Local Data
    TYPE(Matrix_lsr) :: matAT, matBT

    INCLUDE "sparse_includes/SparseBranch.f90"
  END SUBROUTINE SparseBranch_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Calculates the matrix product if using the sparse-sparse algorithm.
  PURE SUBROUTINE SparseBranch_lsc(matA, matB, matC, IsATransposed, &
       & IsBTransposed, alpha, threshold, blocked_memory_pool)
    !> Matrix A.
    TYPE(Matrix_lsc), INTENT(IN)  :: matA
    !> Matrix B
    TYPE(Matrix_lsc), INTENT(IN)  :: matB
    !> matC = alpha*matA*op( matB ) + beta*matC.
    TYPE(Matrix_lsc), INTENT(INOUT) :: matC
    !> True if A is transposed.
    LOGICAL, INTENT(IN) :: IsATransposed
    !> True if B is transposed.
    LOGICAL, INTENT(IN) :: IsBTransposed
    !> Scaling value.
    REAL(NTREAL), INTENT(IN) :: alpha
    !> Threshold for flushing values.
    REAL(NTREAL), INTENT(IN) :: threshold
    !> Memory pool.
    TYPE(MatrixMemoryPool_lc), INTENT(INOUT) :: blocked_memory_pool
    !! Local Data
    TYPE(Matrix_lsc) :: matAT, matBT

    INCLUDE "sparse_includes/SparseBranch.f90"
  END SUBROUTINE SparseBranch_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Calculate the matrix product using the dense-dense algorithm.
  SUBROUTINE DenseBranch_lsr(matA, matB, matC, IsATransposed, IsBTransposed, &
       & alpha, threshold)
    !> Matrix A.
    TYPE(Matrix_lsr), INTENT(IN)  :: matA
    !> Matrix B
    TYPE(Matrix_lsr), INTENT(IN)  :: matB
    !> matC = alpha*matA*op( matB ) + beta*matC.
    TYPE(Matrix_lsr), INTENT(INOUT) :: matC
    !> True if A is transposed.
    LOGICAL, INTENT(IN) :: IsATransposed
    !> True if B is transposed.
    LOGICAL, INTENT(IN) :: IsBTransposed
    !> Scaling value.
    REAL(NTREAL), INTENT(IN) :: alpha
    !> Threshold for flushing values.
    REAL(NTREAL), INTENT(IN) :: threshold
    !! Local Data
    TYPE(Matrix_ldr) :: untransposedMatA
    TYPE(Matrix_ldr) :: untransposedMatB
    TYPE(Matrix_ldr) :: DenseA
    TYPE(Matrix_ldr) :: DenseB
    TYPE(Matrix_ldr) :: DenseC

    INCLUDE "sparse_includes/DenseBranch.f90"
  END SUBROUTINE DenseBranch_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Calculate the matrix product using the dense-dense algorithm.
  SUBROUTINE DenseBranch_lsc(matA, matB, matC, IsATransposed, IsBTransposed, &
       & alpha, threshold)
    !> Matrix A.
    TYPE(Matrix_lsc), INTENT(IN)  :: matA
    !> Matrix B
    TYPE(Matrix_lsc), INTENT(IN)  :: matB
    !> matC = alpha*matA*op( matB ) + beta*matC.
    TYPE(Matrix_lsc), INTENT(INOUT) :: matC
    !> True if A is transposed.
    LOGICAL, INTENT(IN) :: IsATransposed
    !> True if B is transposed.
    LOGICAL, INTENT(IN) :: IsBTransposed
    !> Scaling value.
    REAL(NTREAL), INTENT(IN) :: alpha
    !> Threshold for flushing values.
    REAL(NTREAL), INTENT(IN) :: threshold
    !! Local Data
    TYPE(Matrix_ldc) :: untransposedMatA
    TYPE(Matrix_ldc) :: untransposedMatB
    TYPE(Matrix_ldc) :: DenseA
    TYPE(Matrix_ldc) :: DenseB
    TYPE(Matrix_ldc) :: DenseC

    INCLUDE "sparse_includes/DenseBranch.f90"
  END SUBROUTINE DenseBranch_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Multiplies a single block fo sparse-sparse.
  PURE SUBROUTINE MultiplyBlock_lsr(matAT,matBT,memorypool)
    !> Matrix A, already transposed.
    TYPE(Matrix_lsr), INTENT(IN)  :: matAT
    !> Matrix B, already transposed.
    TYPE(Matrix_lsr), INTENT(IN)  :: matBT
    !> Memory pool to multiply into.
    TYPE(MatrixMemoryPool_lr), INTENT(INOUT) :: memorypool
    !! Temp Variables
    REAL(NTREAL) :: temp_value_a, temp_value_b, temp_value_c

    INCLUDE "sparse_includes/MultiplyBlock.f90"
  END SUBROUTINE MultiplyBlock_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Multiplies a single block fo sparse-sparse.
  PURE SUBROUTINE MultiplyBlock_lsc(matAT,matBT,memorypool)
    !> Matrix A, already transposed.
    TYPE(Matrix_lsc), INTENT(IN)  :: matAT
    !> Matrix B, already transposed.
    TYPE(Matrix_lsc), INTENT(IN)  :: matBT
    !> Memory pool to multiply into.
    TYPE(MatrixMemoryPool_lc), INTENT(INOUT) :: memorypool
    !! Temp Variables
    COMPLEX(NTCOMPLEX) :: temp_value_a, temp_value_b, temp_value_c

    INCLUDE "sparse_includes/MultiplyBlock.f90"
  END SUBROUTINE MultiplyBlock_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Prunes out the values of the hash table into the matrix.
  PURE SUBROUTINE PruneList_lsr(memorypool,alpha,threshold, mat_c_columns, &
       & mat_c_rows, matAB)
    !> Memory pool to prune from.
    TYPE(MatrixMemoryPool_lr), INTENT(INOUT) :: memorypool
    !> Scaling value.
    REAL(NTREAL), INTENT(IN) :: alpha
    !> Threshold for flushing values to zero.
    REAL(NTREAL), INTENT(IN) :: threshold
    !> Size of the matrix we computed (columns).
    INTEGER, INTENT(IN) :: mat_c_columns
    !> Size of the matrix we computed (rows).
    INTEGER, INTENT(IN) :: mat_c_rows
    !> Sparse matrix to prune out into.
    TYPE(Matrix_lsr), INTENT(INOUT) :: matAB
    !! Local data
    REAL(NTREAL) :: working_value
    TYPE(TripletList_r) :: unsorted_pruned_list
    TYPE(TripletList_r) :: sorted_pruned_list

    INCLUDE "sparse_includes/PruneList.f90"
  END SUBROUTINE PruneList_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Prunes out the values of the hash table into the matrix.
  PURE SUBROUTINE PruneList_lsc(memorypool,alpha,threshold, &
       & mat_c_columns, mat_c_rows, matAB)
    !> Memory pool to prune from.
    TYPE(MatrixMemoryPool_lc), INTENT(INOUT) :: memorypool
    !> Scaling value.
    REAL(NTREAL), INTENT(IN) :: alpha
    !> Threshold for flushing values to zero.
    REAL(NTREAL), INTENT(IN) :: threshold
    !> Size of the matrix we computed (columns).
    INTEGER, INTENT(IN) :: mat_c_columns
    !> Size of the matrix we computed (rows).
    INTEGER, INTENT(IN) :: mat_c_rows
    !> Sparse matrix to prune out into.
    TYPE(Matrix_lsc), INTENT(INOUT) :: matAB
    !! Local data
    COMPLEX(NTCOMPLEX) :: working_value
    TYPE(TripletList_c) :: unsorted_pruned_list
    TYPE(TripletList_c) :: sorted_pruned_list

    INCLUDE "sparse_includes/PruneList.f90"
  END SUBROUTINE PruneList_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE SMatrixAlgebraModule
