!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for performing linear algebra using sparse matrices.
MODULE SparseMatrixAlgebraModule
  USE DataTypesModule, ONLY : NTREAL
  USE DenseMatrixModule, ONLY : ConstructDenseFromSparse, &
       & ConstructSparseFromDense, MultiplyDense
  USE MatrixMemoryPoolModule, ONLY : MatrixMemoryPool_t, &
       & ConstructMatrixMemoryPool, DestructMatrixMemoryPool, &
       & CheckMemoryPoolValidity, SetPoolSparsity
  USE SparseMatrixModule, ONLY: SparseMatrix_t, ConstructEmptySparseMatrix, &
       & DestructSparseMatrix, ConstructFromTripletList, CopySparseMatrix, &
       & TransposeSparseMatrix, PrintSparseMatrix
  USE SparseVectorModule, ONLY : AddSparseVectors, DotSparseVectors, &
       & PairwiseMultiplyVectors
  USE TimerModule, ONLY : StartTimer, StopTimer
  USE TripletListModule, ONLY: TripletList_t, SortTripletList, &
       & ConstructTripletList, DestructTripletList
  USE TripletModule, ONLY : Triplet_t
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Linear Algebra
  PUBLIC :: ScaleSparseMatrix
  PUBLIC :: IncrementSparseMatrix
  PUBLIC :: DotSparseMatrix
  PUBLIC :: PairwiseMultiplySparseMatrix
  PUBLIC :: Gemm
  PUBLIC :: SparseMatrixNorm
  PUBLIC :: SparseMatrixGrandSum
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a sparse matrix by a constant.
  !! @param[inout] matA Matrix A.
  !! @param[in] constant scale factor.
  PURE SUBROUTINE ScaleSparseMatrix(matA,constant)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(INOUT) :: matA
    REAL(NTREAL), INTENT(IN) :: constant

    matA%values = constant * matA%values
  END SUBROUTINE ScaleSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Matrix B = alpha*Matrix A + Matrix B (AXPY).
  !! This will utilize the sparse vector addition routine.
  !! @param[in] matA Matrix A.
  !! @param[in,out] matB Matrix B.
  !! @param[in] alpha_in multiplier (optional, default=1.0)
  !! @param[in] threshold_in for flushing values to zero. (Optional, default=0).
  PURE SUBROUTINE IncrementSparseMatrix(matA, matB, alpha_in, threshold_in)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(IN)  :: matA
    TYPE(SparseMatrix_t), INTENT(INOUT) :: matB
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    !! Counter Variables
    INTEGER :: outer_counter
    INTEGER :: elements_per_inner_a, elements_per_inner_b
    INTEGER :: total_counter_a, total_counter_b, total_counter_c
    !! Temporary Variables
    TYPE(SparseMatrix_t) :: matC
    INTEGER :: indices_added_into_c
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: threshold
    INTEGER :: size_of_a, size_of_b

    !! Process Optional Parameters
    IF (.NOT. PRESENT(alpha_in)) THEN
       alpha = 1.0d+0
    ELSE
       alpha = alpha_in
    END IF
    IF (.NOT. PRESENT(threshold_in)) THEN
       threshold = 0.0d+0
    ELSE
       threshold = threshold_in
    END IF

    size_of_a = matA%outer_index(matA%columns+1)

    !! Allocate sufficient space for matC
    CALL ConstructEmptySparseMatrix(matC, matA%columns, matA%rows)
    IF (ALLOCATED(matB%values)) THEN
       size_of_b = matB%outer_index(matB%columns+1)
       ALLOCATE(matC%inner_index(size_of_a+size_of_b))
       ALLOCATE(matC%values(size_of_a+size_of_b))
    ELSE
       ALLOCATE(matC%inner_index(size_of_a))
       ALLOCATE(matC%values(size_of_a))
    END IF

    !! Perform loops
    total_counter_a = 1
    total_counter_b = 1
    total_counter_c = 1
    DO outer_counter = 1, matA%columns
       !! Inner counters
       elements_per_inner_a = matA%outer_index(outer_counter+1) - &
            & matA%outer_index(outer_counter)
       elements_per_inner_b = matB%outer_index(outer_counter+1) - &
            & matB%outer_index(outer_counter)
       CALL AddSparseVectors(&
            matA%inner_index(total_counter_a:total_counter_a+elements_per_inner_a-1),&
            matA%values(total_counter_a:total_counter_a+elements_per_inner_a-1),&
            matB%inner_index(total_counter_b:total_counter_b+elements_per_inner_b-1),&
            matB%values(total_counter_b:total_counter_b+elements_per_inner_b-1),&
            matC%inner_index(total_counter_c:),matC%values(total_counter_c:),&
            indices_added_into_c, alpha, threshold)
       matC%outer_index(outer_counter+1) = matC%outer_index(outer_counter)+&
            & indices_added_into_c
       total_counter_a = total_counter_a + elements_per_inner_a
       total_counter_b = total_counter_b + elements_per_inner_b
       total_counter_c = total_counter_c + indices_added_into_c
    END DO

    !! Cleanup
    CALL DestructSparseMatrix(matB)
    CALL ConstructEmptySparseMatrix(matB,matC%columns,matC%rows)
    matB%outer_index = matC%outer_index
    ALLOCATE(matB%inner_index(matC%outer_index(matC%columns+1)))
    ALLOCATE(matB%values(matC%outer_index(matC%columns+1)))
    matB%inner_index = matC%inner_index(:matC%outer_index(matC%columns+1))
    matB%values = matC%values(:matC%outer_index(matC%columns+1))
    CALL DestructSparseMatrix(matC)
  END SUBROUTINE IncrementSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Pairwise Multiply two matrices.
  !! This will utilize the sparse vector pairwise routine.
  !! @param[in] matA Matrix A.
  !! @param[in] matB Matrix B.
  !! @param[in,out] matC = MatA mult MatB.
  PURE SUBROUTINE PairwiseMultiplySparseMatrix(matA, matB, matC)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(IN)  :: matA
    TYPE(SparseMatrix_t), INTENT(IN) :: matB
    TYPE(SparseMatrix_t), INTENT(INOUT) :: matC
    !! Counter Variables
    INTEGER :: outer_counter
    INTEGER :: elements_per_inner_a, elements_per_inner_b
    INTEGER :: total_counter_a, total_counter_b, total_counter_c
    !! Temporary Variables
    TYPE(SparseMatrix_t) :: TempMat
    INTEGER :: indices_added_into_c
    INTEGER :: size_of_a, size_of_b

    CALL ConstructEmptySparseMatrix(TempMat, matA%columns, matA%rows)
    size_of_a = matA%outer_index(matA%columns+1)
    size_of_b = matB%outer_index(matB%columns+1)
    ALLOCATE(TempMat%inner_index(MIN(size_of_a,size_of_b)))
    ALLOCATE(TempMat%values(MIN(size_of_a,size_of_b)))

    !! Perform loops
    total_counter_a = 1
    total_counter_b = 1
    total_counter_c = 1
    DO outer_counter = 1, matA%columns
       !! Inner counters
       elements_per_inner_a = matA%outer_index(outer_counter+1) - &
            & matA%outer_index(outer_counter)
       elements_per_inner_b = matB%outer_index(outer_counter+1) - &
            & matB%outer_index(outer_counter)
       CALL PairwiseMultiplyVectors(&
            matA%inner_index(total_counter_a:total_counter_a+elements_per_inner_a-1),&
            matA%values(total_counter_a:total_counter_a+elements_per_inner_a-1),&
            matB%inner_index(total_counter_b:total_counter_b+elements_per_inner_b-1),&
            matB%values(total_counter_b:total_counter_b+elements_per_inner_b-1),&
            TempMat%inner_index(total_counter_c:),TempMat%values(total_counter_c:),&
            indices_added_into_c)
       TempMat%outer_index(outer_counter+1) = TempMat%outer_index(outer_counter)+&
            & indices_added_into_c
       total_counter_a = total_counter_a + elements_per_inner_a
       total_counter_b = total_counter_b + elements_per_inner_b
       total_counter_c = total_counter_c + indices_added_into_c
    END DO

    !! Cleanup
    CALL DestructSparseMatrix(matC)
    CALL ConstructEmptySparseMatrix(matC,TempMat%columns,TempMat%rows)
    matC%outer_index = TempMat%outer_index
    ALLOCATE(matC%inner_index(TempMat%outer_index(TempMat%columns+1)))
    ALLOCATE(matC%values(TempMat%outer_index(TempMat%columns+1)))
    matC%inner_index = TempMat%inner_index(:TempMat%outer_index(TempMat%columns+1))
    matC%values = TempMat%values(:TempMat%outer_index(TempMat%columns+1))
    CALL DestructSparseMatrix(TempMat)
  END SUBROUTINE PairwiseMultiplySparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Product = sum(MatA[ij]*MatB[ij])
  !! @param[in] matA Matrix A.
  !! @param[in] matB Matrix B.
  !! @result product
  PURE FUNCTION DotSparseMatrix(matA, matB) RESULT(product)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(IN) :: matA
    TYPE(SparseMatrix_t), INTENT(IN) :: matB
    REAL(NTREAL) :: product
    !! Local Variables
    TYPE(SparseMatrix_t) :: matC

    CALL PairwiseMultiplySparseMatrix(matA,matB,matC)

    product = SparseMatrixGrandSum(matC)
    CALL DestructSparseMatrix(matC)
  END FUNCTION DotSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Multiply two matrices together, and add to the third.
  !! C := alpha*matA*op( matB ) + beta*matC
  !! @param[in] matA Matrix A.
  !! @param[in] matB Matrix B.
  !! @param[out] matC = alpha*matA*op( matB ) + beta*matC.
  !! @param[in] IsATransposed_in true if A is already transposed.
  !! @param[in] IsBTransposed_in true if B is already transposed.
  !! @param[in] alpha_in scales the multiplication.
  !! @param[in] beta_in scales matrix we sum on to.
  !! @param[in] threshold_in for flushing values to zero. Default value is 0.0.
  !! @param[inout] blocked_memory_pool_in an optional memory pool for doing the
  !! calculation.
  SUBROUTINE Gemm(matA, matB, matC, IsATransposed_in, IsBTransposed_in, &
       & alpha_in, beta_in, threshold_in, blocked_memory_pool_in)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(IN)  :: matA
    TYPE(SparseMatrix_t), INTENT(IN)  :: matB
    TYPE(SparseMatrix_t), INTENT(INOUT) :: matC
    LOGICAL, OPTIONAL, INTENT(IN) :: IsATransposed_in
    LOGICAL, OPTIONAL, INTENT(IN) :: IsBTransposed_in
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: beta_in
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    TYPE(MatrixMemoryPool_t), OPTIONAL, &
         & INTENT(INOUT), TARGET :: blocked_memory_pool_in
    !! Intermediate Data
    TYPE(SparseMatrix_t) :: matAB
    LOGICAL :: IsATransposed, IsBTransposed
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: beta
    REAL(NTREAL) :: threshold
    TYPE(MatrixMemoryPool_t) :: blocked_memory_pool
    !! Counters and temporary data
    INTEGER :: mat_c_columns, mat_c_rows
    !! For Efficiency Purposes
    REAL(NTREAL) :: sparsity_a, sparsity_b
    REAL(NTREAL) :: sparsity_estimate
    LOGICAL :: pool_flag

    !! Process Optional Parameters
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
    IF (.NOT. PRESENT(IsATransposed_in)) THEN
       IsATransposed = .FALSE.
    ELSE
       IsATransposed = IsATransposed_in
    END IF
    IF (.NOT. PRESENT(IsBTransposed_in)) THEN
       IsBTransposed = .FALSE.
    ELSE
       IsBTransposed = IsBTransposed_in
    END IF
    IF (.NOT. PRESENT(threshold_in)) THEN
       threshold = 0.0
    ELSE
       threshold = threshold_in
    END IF

    !! Storage details for result matrix
    IF (IsATransposed) THEN
       mat_c_rows = matA%columns
    ELSE
       mat_c_rows = matA%rows
    END IF
    IF (IsBTransposed) THEN
       mat_c_columns = matB%rows
    ELSE
       mat_c_columns = matB%columns
    END IF

    !! Initialization of Memory
    sparsity_a = DBLE(SIZE(matA%values))/(matA%rows*matA%columns)
    sparsity_b = DBLE(SIZE(matB%values))/(matB%rows*matB%columns)
    sparsity_estimate = 4*MAX(sparsity_a,sparsity_b)
    IF (sparsity_estimate > 1.0) THEN
       sparsity_estimate = 1.0
    ELSE IF (sparsity_estimate < 1e-8) THEN
       sparsity_estimate = 1e-8
    END IF

    !! Decide whether to do dense or sparse version.
    IF (MIN(sparsity_a, sparsity_b) .GT. 0.3) THEN
       CALL DenseBranch(matA, matB, matAB, IsATransposed, IsBTransposed, &
            & alpha, threshold)
    ELSE
       !! Setup the memory pool
       IF (.NOT. PRESENT(blocked_memory_pool_in)) THEN
          CALL ConstructMatrixMemoryPool(blocked_memory_pool,mat_c_columns, &
               & mat_c_rows, sparsity_estimate)
          pool_flag = .FALSE.
       ELSEIF (.NOT. CheckMemoryPoolValidity(blocked_memory_pool_in, &
            & mat_c_columns, mat_c_rows)) THEN
          CALL DestructMatrixMemoryPool(blocked_memory_pool_in)
          CALL ConstructMatrixMemoryPool(blocked_memory_pool_in,mat_c_columns, &
               & mat_c_rows, sparsity_estimate)
          pool_flag = .TRUE.
       ELSE
          CALL SetPoolSparsity(blocked_memory_pool_in, sparsity_estimate)
          pool_flag = .TRUE.
       END IF
       !! Multiply
       IF (pool_flag) THEN
          CALL SparseBranch(matA, matB, matAB, IsATransposed, IsBTransposed, &
               & alpha, threshold, blocked_memory_pool_in)
       ELSE
          CALL SparseBranch(matA, matB, matAB, IsATransposed, IsBTransposed, &
               & alpha, threshold, blocked_memory_pool)
       END IF
    END IF

    !! Handle the add part of GEMM
    IF (PRESENT(beta_in)) THEN
       IF (ABS(beta_in) .GT. 0) THEN
          CALL ScaleSparseMatrix(matC,beta)
          CALL IncrementSparseMatrix(matAB,matC)
       ELSE
          CALL CopySparseMatrix(matAB,matC)
       END IF
    ELSE
       CALL CopySparseMatrix(matAB,matC)
    END IF

    CALL DestructSparseMatrix(matAB)
    CALL DestructMatrixMemoryPool(blocked_memory_pool)
  END SUBROUTINE Gemm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the norm of a sparse matrix along the columns.
  !! @param[in] this the matrix to compute the norm of.
  !! @param[out] norm_per_column the norm value for each column in this matrix.
  PURE SUBROUTINE SparseMatrixNorm(this, norm_per_column)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(IN) :: this
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: norm_per_column
    !! Local Data
    INTEGER :: outer_counter, inner_counter
    INTEGER :: elements_per_inner
    REAL(NTREAL) :: temp_value

    !! Allocate Space For Result
    ALLOCATE(norm_per_column(this%columns))
    norm_per_column = 0

    !! Iterate Over Local Data
    DO outer_counter = 1, this%columns
       elements_per_inner = this%outer_index(outer_counter+1) - &
            & this%outer_index(outer_counter)
       DO inner_counter = 1, elements_per_inner
          temp_value = this%values(this%outer_index(outer_counter)+ &
               & inner_counter)
          norm_per_column(outer_counter) = norm_per_column(outer_counter) + &
               & ABS(temp_value)
       END DO
    END DO
  END SUBROUTINE SparseMatrixNorm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sum the elements of a matrix
  !! @param[in] this the matrix to sum
  !! @result sum_value the sum of the matrix elements
  PURE FUNCTION SparseMatrixGrandSum(this) RESULT(sum_value)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(IN) :: this
    REAL(NTREAL) :: sum_value

    sum_value = SUM(this%values)

  END FUNCTION SparseMatrixGrandSum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE SUBROUTINE SparseBranch(matA, matB, matC, IsATransposed, IsBTransposed, &
       & alpha, threshold, blocked_memory_pool)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(IN)  :: matA
    TYPE(SparseMatrix_t), INTENT(IN)  :: matB
    TYPE(SparseMatrix_t), INTENT(INOUT) :: matC
    LOGICAL, INTENT(IN) :: IsATransposed
    LOGICAL, INTENT(IN) :: IsBTransposed
    REAL(NTREAL), INTENT(IN) :: alpha
    REAL(NTREAL), INTENT(IN) :: threshold
    TYPE(MatrixMemoryPool_t), INTENT(INOUT) :: blocked_memory_pool
    !! Local Data
    TYPE(SparseMatrix_t) :: matAT, matBT

    !! Block A and B
    IF (.NOT. IsATransposed) THEN
       CALL TransposeSparseMatrix(matA,matAT)
    END IF
    IF (.NOT. IsBTransposed) THEN
       CALL TransposeSparseMatrix(matB,matBT)
    END IF

    IF (IsATransposed .AND. IsBTransposed) THEN
       CALL MultiplyBlock(matA, matB, blocked_memory_pool)
    ELSEIF (IsATransposed) THEN
       CALL MultiplyBlock(matA, matBT, blocked_memory_pool)
    ELSEIF (IsBTransposed) THEN
       CALL MultiplyBlock(matAT, matB, blocked_memory_pool)
    ELSE
       CALL MultiplyBlock(matAT, matBT, blocked_memory_pool)
    END IF

    !! Go from triplets to return matrix
    CALL PruneList(blocked_memory_pool, alpha, threshold, &
         & blocked_memory_pool%columns, blocked_memory_pool%rows, matC)
  END SUBROUTINE SparseBranch
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE SUBROUTINE DenseBranch(matA, matB, matC, IsATransposed, IsBTransposed, &
       & alpha, threshold)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(IN)  :: matA
    TYPE(SparseMatrix_t), INTENT(IN)  :: matB
    TYPE(SparseMatrix_t), INTENT(INOUT) :: matC
    LOGICAL, INTENT(IN) :: IsATransposed
    LOGICAL, INTENT(IN) :: IsBTransposed
    REAL(NTREAL), INTENT(IN) :: alpha
    REAL(NTREAL), INTENT(IN) :: threshold
    !! Local Data
    TYPE(SparseMatrix_t) :: untransposedMatA
    TYPE(SparseMatrix_t) :: untransposedMatB
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: DenseA
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: DenseB
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: DenseC

    !! Handle Transposed Case
    IF (IsATransposed) THEN
       CALL TransposeSparseMatrix(matA,untransposedMatA)
    ELSE
       CALL CopySparseMatrix(matA,untransposedMatA)
    END IF
    IF (IsBTransposed) THEN
       CALL TransposeSparseMatrix(matB,untransposedMatB)
    ELSE
       CALL CopySparseMatrix(matB,untransposedMatB)
    END IF

    !! Convert Forward
    ALLOCATE(DenseA(untransposedMatA%rows,untransposedMatA%columns))
    ALLOCATE(DenseB(untransposedMatB%rows,untransposedMatB%columns))
    ALLOCATE(DenseC(untransposedMatA%rows,untransposedMatB%columns))
    CALL ConstructDenseFromSparse(untransposedMatA, DenseA)
    CALL ConstructDenseFromSparse(untransposedMatB, DenseB)

    !! Multiply
    CALL MultiplyDense(DenseA, DenseB, DenseC)

    !! Convert Back
    CALL ConstructSparseFromDense(DenseC, matC, threshold)
    CALL ScaleSparseMatrix(matC,alpha)

    DEALLOCATE(DenseA)
    DEALLOCATE(DenseB)
    DEALLOCATE(DenseC)
  END SUBROUTINE DenseBranch
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE SUBROUTINE MultiplyBlock(matAT,matBT,memorypool)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(IN)  :: matAT
    TYPE(SparseMatrix_t), INTENT(IN)  :: matBT
    TYPE(MatrixMemoryPool_t), INTENT(INOUT) :: memorypool
    !! Temp Variables
    REAL(NTREAL) :: temp_value_a, temp_value_b, temp_value_c
    INTEGER :: temp_inserted_values
    INTEGER :: temp_index_a, temp_index_b
    INTEGER :: elements_per_inner_a
    INTEGER :: elements_per_inner_b
    LOGICAL :: is_set
    !! Counters
    INTEGER :: outer_counter, inner_counter_a, inner_counter_b

    !! Multiply
    DO outer_counter = 1, matAT%columns
       elements_per_inner_a = matAT%outer_index(outer_counter+1) - &
            & matAT%outer_index(outer_counter)
       DO inner_counter_a = 1, elements_per_inner_a
          temp_value_a = matAT%values(matAT%outer_index(outer_counter)+ &
               & inner_counter_a)
          temp_index_a = matAT%inner_index(matAT%outer_index(outer_counter)+ &
               & inner_counter_a)
          elements_per_inner_b = matBT%outer_index(temp_index_a+1) - &
               & matBT%outer_index(temp_index_a)
          DO inner_counter_b = 1, elements_per_inner_b
             temp_index_b = matBT%inner_index(matBT%outer_index(temp_index_a)+ &
                  & inner_counter_b)
             temp_value_b = matBT%values(matBT%outer_index(temp_index_a)+ &
                  & inner_counter_b)
             temp_value_c = memorypool%value_array(temp_index_b,outer_counter)
             is_set = memorypool%dirty_array(temp_index_b,outer_counter)
             IF (is_set .EQV. .FALSE.) THEN
                memorypool%dirty_array(temp_index_b,outer_counter) = .TRUE.
                temp_inserted_values = memorypool%inserted_per_bucket(&
                     & (temp_index_b-1)/memorypool%hash_size+1,outer_counter) + 1
                memorypool%inserted_per_bucket(&
                     & (temp_index_b-1)/memorypool%hash_size+1,outer_counter) = &
                     & temp_inserted_values
                memorypool%hash_index(temp_inserted_values+ &
                     & ((temp_index_b-1)/memorypool%hash_size)*memorypool%hash_size, &
                     & outer_counter) = temp_index_b
             END IF
             memorypool%value_array(temp_index_b,outer_counter) = &
                  & temp_value_c + temp_value_a*temp_value_b
          END DO
       END DO
    END DO
  END SUBROUTINE MultiplyBlock
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE SUBROUTINE PruneList(memorypool,alpha,threshold, &
       & mat_c_columns, mat_c_rows, matAB)
    !! Parameters
    TYPE(MatrixMemoryPool_t), INTENT(INOUT) :: memorypool
    REAL(NTREAL), INTENT(IN) :: alpha
    REAL(NTREAL), INTENT(IN) :: threshold
    INTEGER, INTENT(IN) :: mat_c_columns
    INTEGER, INTENT(IN) :: mat_c_rows
    TYPE(SparseMatrix_t), INTENT(INOUT) :: matAB
    !! Local data
    INTEGER :: row_counter_c, column_counter_c, hash_counter
    REAL(NTREAL) :: working_value
    INTEGER :: working_column
    INTEGER :: temp_values_per_hash
    TYPE(TripletList_t) :: unsorted_pruned_list
    TYPE(TripletList_t) :: sorted_pruned_list
    INTEGER :: pruned_counter

    pruned_counter = 1
    DO row_counter_c = 1, mat_c_rows
       DO column_counter_c = 1, (mat_c_columns-1)/memorypool%hash_size+1
          !! Sort the elements in a hash
          temp_values_per_hash = memorypool%inserted_per_bucket(&
               & column_counter_c,row_counter_c)
          memorypool%inserted_per_bucket(column_counter_c,row_counter_c) = 0
          !! Copy them
          DO hash_counter=1,temp_values_per_hash
             working_column = memorypool%hash_index(hash_counter+ &
                  & (column_counter_c-1)*memorypool%hash_size, row_counter_c)
             working_value = memorypool%value_array(working_column,row_counter_c)
             memorypool%value_array(working_column,row_counter_c) = 0
             memorypool%dirty_array(working_column,row_counter_c) = .FALSE.
             IF (ABS(alpha*working_value) .GT. threshold) THEN
                memorypool%pruned_list(pruned_counter)%point_value = &
                     & alpha*working_value
                memorypool%pruned_list(pruned_counter)%index_column = &
                     & working_column
                memorypool%pruned_list(pruned_counter)%index_row = &
                     & row_counter_c
                pruned_counter = pruned_counter + 1
             END IF
          END DO
       END DO
    END DO
    CALL ConstructTripletList(unsorted_pruned_list,pruned_counter-1)
    unsorted_pruned_list%data = memorypool%pruned_list(1:pruned_counter-1)
    CALL SortTripletList(unsorted_pruned_list,mat_c_columns,sorted_pruned_list,&
         & bubble_in=.TRUE.)
    CALL ConstructFromTripletList(matAB,sorted_pruned_list,mat_c_rows, &
         & mat_c_columns)
    CALL DestructTripletList(sorted_pruned_list)
  END SUBROUTINE PruneList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE SparseMatrixAlgebraModule
