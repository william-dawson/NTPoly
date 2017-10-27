!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Performing Distributed Sparse Matrix Algebra Operations.
MODULE DistributedSparseMatrixAlgebraModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL
  USE DistributedMatrixMemoryPoolModule, ONLY : &
       & DistributedMatrixMemoryPool_t, ConstructDistributedMatrixMemoryPool, &
       & CheckDistributedMemoryPoolValidity
  USE DistributedSparseMatrixModule
  USE GemmTasksModule
  USE MatrixGatherModule, ONLY : GatherHelper_t, GatherSizes, &
       & GatherAndComposeData, GatherAndComposeCleanup, GatherAndSumData, &
       & GatherAndSumCleanup, TestSizeRequest, TestOuterRequest, &
       & TestInnerRequest, TestDataRequest
  USE ProcessGridModule, ONLY : grid_error,&
       & num_process_slices, my_slice, &
       & within_slice_comm, row_comm, column_comm, between_slice_comm, &
       & blocked_column_comm, blocked_row_comm, blocked_between_slice_comm, &
       & number_of_blocks_rows, number_of_blocks_columns, block_multiplier
  USE SparseMatrixAlgebraModule, ONLY : &
       & DotSparseMatrix, PairwiseMultiplySparseMatrix, &
       & SparseMatrixNorm, ScaleSparseMatrix, IncrementSparseMatrix, Gemm, &
       & SparseMatrixGrandSum
  USE SparseMatrixModule, ONLY : SparseMatrix_t, DestructSparseMatrix, &
       & CopySparseMatrix, &
       & TransposeSparseMatrix, ComposeSparseMatrixColumns, MatrixToTripletList
  USE TimerModule, ONLY : StartTimer, StopTimer
  USE TripletListModule, ONLY : TripletList_t, DestructTripletList
  USE ISO_C_BINDING
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Basic Linear Algebra
  PUBLIC :: ComputeSigma
  PUBLIC :: DistributedGemm
  PUBLIC :: DistributedGrandSum
  PUBLIC :: DistributedPairwiseMultiply
  PUBLIC :: DistributedSparseNorm
  PUBLIC :: DotDistributedSparseMatrix
  PUBLIC :: IncrementDistributedSparseMatrix
  PUBLIC :: ScaleDistributedSparseMatrix
  PUBLIC :: Trace
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute sigma for the inversion method.
  !! @todo describe this better.
  !! @param[in] this the matrix to compute the sigma value of.
  !! @param[out] sigma_value sigma.
  SUBROUTINE ComputeSigma(this,sigma_value)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(in) :: this
    REAL(NTREAL), INTENT(out) :: sigma_value
    !! Local Data
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: &
         & column_sigma_contribution
    !! Counters/Temporary
    INTEGER :: inner_counter, outer_counter
    TYPE(SparseMatrix_t) :: merged_local_data

    !! Merge all the local data
    CALL MergeLocalBlocks(this, merged_local_data)

    ALLOCATE(column_sigma_contribution(merged_local_data%columns))
    column_sigma_contribution = 0
    DO outer_counter = 1, merged_local_data%columns
       DO inner_counter = merged_local_data%outer_index(outer_counter), &
            & merged_local_data%outer_index(outer_counter+1)-1
          column_sigma_contribution(outer_counter) = &
               & column_sigma_contribution(outer_counter) + &
               & ABS(merged_local_data%values(inner_counter+1))
       END DO
    END DO
    CALL MPI_Allreduce(MPI_IN_PLACE,column_sigma_contribution,&
         & merged_local_data%columns,MPINTREAL,MPI_SUM, &
         & column_comm, grid_error)
    CALL MPI_Allreduce(MAXVAL(column_sigma_contribution),sigma_value,1, &
         & MPINTREAL,MPI_MAX,row_comm,grid_error)
    sigma_value = 1.0d+0/(sigma_value**2)

    DEALLOCATE(column_sigma_contribution)
    CALL DestructSparseMatrix(merged_local_data)
  END SUBROUTINE ComputeSigma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Multiply two matrices together, and add to the third.
  !! C := alpha*matA*matB+ beta*matC
  !! @param[in] matA Matrix A
  !! @param[in] matB Matrix B
  !! @param[out] matC = alpha*matA*matB + beta*matC
  !! @param[in] alpha_in scales the multiplication
  !! @param[in] beta_in scales matrix we sum on to
  !! @param[in] threshold_in for flushing values to zero. Default value is 0.0.
  !! @param[inout] memory_pool_in a memory pool that can be used for the
  !! calculation.
  SUBROUTINE DistributedGemm(matA,matB,matC,alpha_in,beta_in,threshold_in, &
       & memory_pool_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(in)    :: matA
    TYPE(DistributedSparseMatrix_t), INTENT(in)    :: matB
    TYPE(DistributedSparseMatrix_t), INTENT(inout) :: matC
    REAL(NTREAL), OPTIONAL, INTENT(in) :: alpha_in
    REAL(NTREAL), OPTIONAL, INTENT(in) :: beta_in
    REAL(NTREAL), OPTIONAL, INTENT(in) :: threshold_in
    TYPE(DistributedMatrixMemoryPool_t), OPTIONAL, INTENT(inout) :: &
         & memory_pool_in
    !! Local Versions of Optional Parameter
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: beta
    REAL(NTREAL) :: threshold
    TYPE(DistributedSparseMatrix_t) :: matAB
    !! Temporary Matrices
    TYPE(SparseMatrix_t), DIMENSION(:,:), ALLOCATABLE :: AdjacentABlocks
    TYPE(SparseMatrix_t), DIMENSION(:), ALLOCATABLE :: LocalRowContribution
    TYPE(SparseMatrix_t), DIMENSION(:), ALLOCATABLE :: GatheredRowContribution
    TYPE(SparseMatrix_t), DIMENSION(:), ALLOCATABLE :: GatheredRowContributionT
    TYPE(SparseMatrix_t), DIMENSION(:,:), ALLOCATABLE :: TransposedBBlocks
    TYPE(SparseMatrix_t), DIMENSION(:), ALLOCATABLE :: LocalColumnContribution
    TYPE(SparseMatrix_t), DIMENSION(:), ALLOCATABLE :: GatheredColumnContribution
    TYPE(SparseMatrix_t), DIMENSION(:,:), ALLOCATABLE :: SliceContribution
    !! Communication Helpers
    TYPE(GatherHelper_t), DIMENSION(:), ALLOCATABLE :: row_helper
    TYPE(GatherHelper_t), DIMENSION(:), ALLOCATABLE :: column_helper
    TYPE(GatherHelper_t), DIMENSION(:,:), ALLOCATABLE :: slice_helper
    !! For Iterating Over Local Blocks
    INTEGER :: row_counter, inner_row_counter
    INTEGER :: column_counter, inner_column_counter
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

    CALL StartTimer("GEMM")

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
    !! The threshhold needs to be smaller if we're doing a sliced version
    !! because you might flush a value that would be kept in the summed version.
    IF (num_process_slices .GT. 1) THEN
       working_threshold = threshold/(num_process_slices*1000)
    ELSE
       working_threshold = threshold
    END IF

    !! Construct The Temporary Matrices
    CALL ConstructEmptyDistributedSparseMatrix(matAB, &
         & matA%actual_matrix_dimension)

    ALLOCATE(AdjacentABlocks(number_of_blocks_columns/num_process_slices, &
         & number_of_blocks_rows))
    ALLOCATE(LocalRowContribution(number_of_blocks_rows))
    ALLOCATE(GatheredRowContribution(number_of_blocks_rows))
    ALLOCATE(GatheredRowContributionT(number_of_blocks_rows))

    ALLOCATE(TransposedBBlocks(number_of_blocks_columns, &
         & number_of_blocks_rows/num_process_slices))
    ALLOCATE(LocalColumnContribution(number_of_blocks_columns))
    ALLOCATE(GatheredColumnContribution(number_of_blocks_columns))
    ALLOCATE(SliceContribution(number_of_blocks_columns, &
         & number_of_blocks_rows))

    !! Helpers
    ALLOCATE(row_helper(number_of_blocks_rows))
    ALLOCATE(column_helper(number_of_blocks_columns))
    ALLOCATE(slice_helper(number_of_blocks_columns, &
         & number_of_blocks_rows))

    !! Construct the task queues
    ALLOCATE(ATasks(number_of_blocks_rows))
    DO row_counter=1,number_of_blocks_rows
       ATasks(row_counter) = LocalGatherA
    END DO
    ALLOCATE(BTasks(number_of_blocks_columns))
    DO column_counter=1,number_of_blocks_columns
       BTasks(column_counter) = LocalGatherB
    END DO
    ALLOCATE(ABTasks(number_of_blocks_columns,number_of_blocks_rows))
    DO row_counter=1,number_of_blocks_rows
       DO column_counter=1,number_of_blocks_columns
          ABTasks(column_counter,row_counter) = AwaitingAB
       END DO
    END DO

    !! Setup A Tasks
    duplicate_start_column = my_slice+1
    duplicate_offset_column = num_process_slices

    !! Setup B Tasks
    duplicate_start_row = my_slice+1
    duplicate_offset_row = num_process_slices

    !! Setup AB Tasks
    IF (PRESENT(memory_pool_in)) THEN
       IF (.NOT. CheckDistributedMemoryPoolValidity(memory_pool_in)) THEN
          CALL ConstructDistributedMatrixMemoryPool(memory_pool_in)
       END IF
    END IF

    !! Run A Tasks
    ATasks_completed = 0
    BTasks_completed = 0
    ABTasks_completed = 0
    !$omp PARALLEL
    !$omp MASTER
    DO WHILE (ATasks_completed .LT. SIZE(ATasks) .OR. &
         & BTasks_completed .LT. SIZE(BTasks) .OR. &
         & ABTasks_completed .LT. SIZE(ABTasks))
       DO row_counter=1,number_of_blocks_rows
          SELECT CASE (ATasks(row_counter))
          CASE(LocalGatherA)
             ATasks(row_counter) = TaskRunningA
             !$omp task default(shared), private(inner_column_counter), &
             !$omp& firstprivate(row_counter)
             !! First Align The Data We're Working With
             DO inner_column_counter=1, &
                  & number_of_blocks_columns/num_process_slices
                CALL CopySparseMatrix(matA%local_data(duplicate_start_column+ &
                     & duplicate_offset_column*(inner_column_counter-1), &
                     & row_counter), &
                     & AdjacentABlocks(inner_column_counter,row_counter))
             END DO
             !! Then Do A Local Gather
             CALL ComposeSparseMatrixColumns(AdjacentABlocks(:,row_counter), &
                  & LocalRowContribution(row_counter))
             ATasks(row_counter) = SendSizeA
             !$omp end task
          CASE(SendSizeA)
             !! Then Start A Global Gather
             CALL GatherSizes(LocalRowContribution(row_counter), &
                  & blocked_row_comm(row_counter), row_helper(row_counter))
             ATasks(row_counter) = ComposeA
          CASE(ComposeA)
             IF (TestSizeRequest(row_helper(row_counter))) THEN
                CALL GatherAndComposeData(LocalRowContribution(row_counter), &
                     & blocked_row_comm(row_counter), &
                     & GatheredRowContribution(row_counter), &
                     & row_helper(row_counter))
                ATasks(row_counter) = WaitOuterA
             END IF
          CASE(WaitOuterA)
             IF (TestOuterRequest(row_helper(row_counter))) THEN
                ATasks(row_counter) = WaitInnerA
             END IF
          CASE(WaitInnerA)
             IF (TestInnerRequest(row_helper(row_counter))) THEN
                ATasks(row_counter) = WaitDataA
             END IF
          CASE(WaitDataA)
             IF (TestDataRequest(row_helper(row_counter))) THEN
                ATasks(row_counter) = AdjustIndicesA
             END IF
          CASE(AdjustIndicesA)
             ATasks(row_counter) = TaskRunningA
             !$omp task default(shared), firstprivate(row_counter)
             CALL GatherAndComposeCleanup(LocalRowContribution(row_counter), &
                  & GatheredRowContribution(row_counter), &
                  & row_helper(row_counter))
             CALL TransposeSparseMatrix(GatheredRowContribution(row_counter), &
                  & GatheredRowContributionT(row_counter))
             ATasks(row_counter) = CleanupA
             !$omp end task
          CASE(CleanupA)
             ATasks(row_counter) = FinishedA
             ATasks_completed = ATasks_completed + 1
          END SELECT
       END DO
       !! B Tasks
       DO column_counter=1,number_of_blocks_columns
          SELECT CASE (BTasks(column_counter))
          CASE(LocalGatherB)
             BTasks(column_counter) = TaskRunningB
             !$omp task default(shared), private(inner_row_counter), &
             !$omp& firstprivate(column_counter)
             !! First Transpose The Data We're Working With
             DO inner_row_counter=1,number_of_blocks_rows/num_process_slices
                CALL TransposeSparseMatrix(matB%local_data(column_counter, &
                     & duplicate_start_row+ &
                     & duplicate_offset_row*(inner_row_counter-1)), &
                     & TransposedBBlocks(column_counter,inner_row_counter))
             END DO
             !! Then Do A Local Gather
             CALL ComposeSparseMatrixColumns(&
                  & TransposedBBlocks(column_counter,:),&
                  & LocalColumnContribution(column_counter))
             BTasks(column_counter) = SendSizeB
             !$omp end task
          CASE(SendSizeB)
             !! Then A Global Gather
             CALL GatherSizes(LocalColumnContribution(column_counter), &
                  & blocked_column_comm(column_counter), &
                  & column_helper(column_counter))
             BTasks(column_counter) = LocalComposeB
          CASE(LocalComposeB)
             IF (TestSizeRequest(column_helper(column_counter))) THEN
                CALL GatherAndComposeData( &
                     & LocalColumnContribution(column_counter),&
                     & blocked_column_comm(column_counter), &
                     & GatheredColumnContribution(column_counter), &
                     & column_helper(column_counter))
                BTasks(column_counter) = WaitOuterB
             END IF
          CASE(WaitOuterB)
             IF (TestOuterRequest(column_helper(column_counter))) THEN
                BTasks(column_counter) = WaitInnerB
             END IF
          CASE(WaitInnerB)
             IF (TestInnerRequest(column_helper(column_counter))) THEN
                BTasks(column_counter) = WaitDataB
             END IF
          CASE(WaitDataB)
             IF (TestDataRequest(column_helper(column_counter))) THEN
                BTasks(column_counter) = AdjustIndicesB
             END IF
          CASE(AdjustIndicesB)
             BTasks(column_counter) = TaskRunningB
             !$omp task default(shared), firstprivate(column_counter)
             CALL GatherAndComposeCleanup( &
                  & LocalColumnContribution(column_counter),&
                  & GatheredColumnContribution(column_counter), &
                  & column_helper(column_counter))
             BTasks(column_counter) = CleanupB
             !$omp end task
          CASE(CleanupB)
             BTasks(column_counter) = FinishedB
             BTasks_completed = BTasks_completed + 1
          END SELECT
       END DO
       !! AB Tasks
       DO row_counter=1,number_of_blocks_rows
          DO column_counter=1,number_of_blocks_columns
             SELECT CASE(ABTasks(column_counter,row_counter))
             CASE (AwaitingAB)
                IF (ATasks(row_counter) .EQ. FinishedA .AND. &
                     & BTasks(column_counter) .EQ. FinishedB) THEN
                   ABTasks(column_counter,row_counter) = GemmAB
                END IF
             CASE (GemmAB)
                ABTasks(column_counter,row_counter) = TaskRunningAB
                !$omp task default(shared), &
                !$omp& firstprivate(row_counter,column_counter)
                IF (PRESENT(memory_pool_in)) THEN
                   CALL Gemm(GatheredRowContributionT(row_counter), &
                        & GatheredColumnContribution(column_counter), &
                        & SliceContribution(column_counter,row_counter), &
                        & IsATransposed_in=.TRUE., IsBTransposed_in=.TRUE., &
                        & alpha_in=alpha, threshold_in=working_threshold, &
                        & blocked_memory_pool_in= &
                        & memory_pool_in%grid(column_counter,row_counter))
                ELSE
                   CALL Gemm(GatheredRowContributionT(row_counter), &
                        & GatheredColumnContribution(column_counter), &
                        & SliceContribution(column_counter,row_counter), &
                        & IsATransposed_in=.TRUE., IsBTransposed_in=.TRUE., &
                        & alpha_in=alpha, threshold_in=working_threshold)
                END IF
                ABTasks(column_counter,row_counter) = SendSizeAB
                !$omp end task
             CASE(SendSizeAB)
                CALL GatherSizes(SliceContribution(column_counter,row_counter), &
                     & blocked_between_slice_comm(column_counter,row_counter), &
                     & slice_helper(column_counter,row_counter))
                ABTasks(column_counter,row_counter) = GatherAndSumAB
             CASE (GatherAndSumAB)
                IF (TestSizeRequest(slice_helper(column_counter,row_counter))) THEN
                   CALL GatherAndSumData( &
                        & SliceContribution(column_counter,row_counter), &
                        & blocked_between_slice_comm(column_counter,row_counter), &
                        & slice_helper(column_counter,row_counter))
                   ABTasks(column_counter,row_counter) = WaitOuterAB
                END IF
             CASE (WaitOuterAB)
                IF (TestOuterRequest(slice_helper(column_counter,row_counter))) THEN
                   ABTasks(column_counter,row_counter) = WaitInnerAB
                END IF
             CASE (WaitInnerAB)
                IF (TestInnerRequest(slice_helper(column_counter,row_counter))) THEN
                   ABTasks(column_counter,row_counter) = WaitDataAB
                END IF
             CASE (WaitDataAB)
                IF (TestDataRequest(slice_helper(column_counter,row_counter))) THEN
                   ABTasks(column_counter,row_counter) = LocalSumAB
                END IF
             CASE(LocalSumAB)
                ABTasks(column_counter,row_counter) = TaskRunningAB
                !$omp task default(shared), &
                !$omp& firstprivate(row_counter,column_counter)
                CALL GatherAndSumCleanup( &
                     & SliceContribution(column_counter,row_counter), &
                     & matAB%local_data(column_counter,row_counter), &
                     & threshold, slice_helper(column_counter,row_counter))
                ABTasks(column_counter,row_counter) = CleanupAB
                !$omp end task
             CASE(CleanupAB)
                ABTasks(column_counter,row_counter) = FinishedAB
                ABTasks_completed = ABTasks_completed + 1
             END SELECT
          END DO
       END DO
    END DO
    !$omp end MASTER
    !$omp end PARALLEL

    !! Copy to output matrix.
    IF (beta .EQ. 0.0) THEN
       CALL CopyDistributedSparseMatrix(matAB,matC)
    ELSE
       CALL ScaleDistributedSparseMatrix(MatC,beta)
       CALL IncrementDistributedSparseMatrix(MatAB,MatC)
    END IF

    !! Cleanup
    CALL DestructDistributedSparseMatrix(matAB)
    DEALLOCATE(row_helper)
    DEALLOCATE(column_helper)
    DEALLOCATE(slice_helper)

    !! Deallocate Buffers From A
    DO row_counter=1,number_of_blocks_rows
       DO inner_column_counter=1,number_of_blocks_columns/num_process_slices
          CALL DestructSparseMatrix(AdjacentABlocks(inner_column_counter, &
               & row_counter))
       END DO
       CALL DestructSparseMatrix(LocalRowContribution(row_counter))
       CALL DestructSparseMatrix(GatheredRowContribution(row_counter))
    END DO
    DEALLOCATE(AdjacentABlocks)
    DEALLOCATE(LocalRowContribution)
    DEALLOCATE(GatheredRowContribution)
    !! Deallocate Buffers From B
    DO column_counter=1,number_of_blocks_columns
       DO inner_row_counter=1,number_of_blocks_rows/num_process_slices
          CALL DestructSparseMatrix(TransposedBBlocks(column_counter, &
               & inner_row_counter))
       END DO
       CALL DestructSparseMatrix(LocalColumnContribution(column_counter))
    END DO
    DEALLOCATE(TransposedBBlocks)
    DEALLOCATE(LocalColumnContribution)
    !! Deallocate Buffers From Multiplying The Block
    DO row_counter=1,number_of_blocks_rows
       CALL DestructSparseMatrix(GatheredRowContributionT(row_counter))
    END DO
    DO column_counter=1,number_of_blocks_columns
       CALL DestructSparseMatrix(GatheredColumnContribution(column_counter))
    END DO
    DEALLOCATE(GatheredRowContributionT)
    DEALLOCATE(GatheredColumnContribution)
    !! Deallocate Buffers From Sum
    DO row_counter=1,number_of_blocks_rows
       DO column_counter=1,number_of_blocks_columns
          CALL DestructSparseMatrix(SliceContribution(column_counter,row_counter))
       END DO
    END DO
    DEALLOCATE(SliceContribution)

    CALL StopTimer("GEMM")
  END SUBROUTINE DistributedGemm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sum up the elements in a matrix.
  !! @param[in] matA Matrix A.
  !! @result sum the sum of all elements.
  FUNCTION DistributedGrandSum(matA) RESULT(sum)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(in)  :: matA
    REAL(NTREAL) :: sum
    !! Local Data
    INTEGER :: row_counter, column_counter
    REAL(NTREAL) :: temp

    sum = 0
    DO row_counter = 1, number_of_blocks_rows
       DO column_counter = 1, number_of_blocks_columns
          temp = SparseMatrixGrandSum(matA%local_data(column_counter,row_counter))
          sum = sum + temp
       END DO
    END DO

    !! Sum Among Process Slice
    CALL MPI_Allreduce(MPI_IN_PLACE, sum, 1, MPINTREAL, &
         & MPI_SUM, within_slice_comm, grid_error)
  END FUNCTION DistributedGrandSum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Elementwise multiplication. C_ij = A_ij * B_ij.
  !! @param[in] matA Matrix A.
  !! @param[in] matB Matrix B.
  !! @param[in,out] matC = MatA mult MatB.
  SUBROUTINE DistributedPairwiseMultiply(matA, matB, matC)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(in)  :: matA
    TYPE(DistributedSparseMatrix_t), INTENT(in)  :: matB
    TYPE(DistributedSparseMatrix_t), INTENT(inout)  :: matC
    !! Local Data
    INTEGER :: row_counter, column_counter

    CALL ConstructEmptyDistributedSparseMatrix(matC, &
         & matA%actual_matrix_dimension)

    !$omp parallel
    !$omp do collapse(2)
    DO row_counter = 1, number_of_blocks_rows
       DO column_counter = 1, number_of_blocks_columns
          CALL PairwiseMultiplySparseMatrix( &
               & matA%local_data(column_counter,row_counter), &
               & matB%local_data(column_counter,row_counter), &
               & matC%local_data(column_counter,row_counter))
       END DO
    END DO
    !$omp end do
    !$omp end parallel
  END SUBROUTINE DistributedPairwiseMultiply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the norm of a distributed sparse matrix along the rows.
  !! @param[in] this the matrix to compute the norm of.
  !! @return the norm value of the full distributed sparse matrix.
  FUNCTION DistributedSparseNorm(this) RESULT(norm_value)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(in) :: this
    REAL(NTREAL) :: norm_value
    !! Local Data
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: local_norm
    TYPE(SparseMatrix_t) :: merged_local_data

    !! Merge all the local data
    CALL MergeLocalBlocks(this, merged_local_data)

    !! Sum Along Columns
    CALL SparseMatrixNorm(merged_local_data,local_norm)
    CALL MPI_Allreduce(MPI_IN_PLACE,local_norm,SIZE(local_norm), &
         & MPINTREAL, MPI_SUM,column_comm,grid_error)

    !! Find Max Value Amonst Columns
    norm_value = MAXVAL(local_norm)
    CALL MPI_Allreduce(MPI_IN_PLACE,norm_value,1,MPINTREAL,MPI_MAX, &
         & row_comm, grid_error)

    CALL DestructSparseMatrix(merged_local_data)
    DEALLOCATE(local_norm)
  END FUNCTION DistributedSparseNorm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> product = dot(Matrix A,Matrix B)
  !! @param[in] matA Matrix A.
  !! @param[in,out] matB Matrix B.
  !! @result product the dot product.
  FUNCTION DotDistributedSparseMatrix(matA, matB) RESULT(product)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(in)  :: matA
    TYPE(DistributedSparseMatrix_t), INTENT(in)  :: matB
    REAL(NTREAL) :: product
    !! Local Data
    TYPE(DistributedSparseMatrix_t)  :: matC

    CALL DistributedPairwiseMultiply(matA,matB,matC)
    product = DistributedGrandSum(matC)
    CALL DestructDistributedSparseMatrix(matC)
  END FUNCTION DotDistributedSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Matrix B = alpha*Matrix A + Matrix B (AXPY)
  !! This will utilize the sparse vector increment routine.
  !! @param[in] matA Matrix A.
  !! @param[in,out] matB Matrix B.
  !! @param[in] alpha_in multiplier. Default value is 1.0
  !! @param[in] threshold_in for flushing values to zero. Default value is 0.0
  SUBROUTINE IncrementDistributedSparseMatrix(matA, matB, alpha_in,threshold_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(in)  :: matA
    TYPE(DistributedSparseMatrix_t), INTENT(inout)  :: matB
    REAL(NTREAL), OPTIONAL, INTENT(in) :: alpha_in
    REAL(NTREAL), OPTIONAL, INTENT(in) :: threshold_in
    !! Local Data
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: threshold
    INTEGER :: row_counter, column_counter

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

    !$omp parallel
    !$omp do collapse(2)
    DO row_counter = 1, number_of_blocks_rows
       DO column_counter = 1, number_of_blocks_columns
          CALL IncrementSparseMatrix(matA%local_data(column_counter,row_counter),&
               & matB%local_data(column_counter,row_counter), alpha, threshold)
       END DO
    END DO
    !$omp end do
    !$omp end parallel

  END SUBROUTINE IncrementDistributedSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a distributed sparse matrix by a constant.
  !! @param[inout] this Matrix to scale.
  !! @param[in] constant scale factor.
  SUBROUTINE ScaleDistributedSparseMatrix(this,constant)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(inout) :: this
    REAL(NTREAL), INTENT(in) :: constant
    INTEGER :: row_counter, column_counter

    !$omp parallel
    !$omp do collapse(2)
    DO row_counter = 1, number_of_blocks_rows
       DO column_counter = 1, number_of_blocks_columns
          CALL ScaleSparseMatrix(this%local_data(column_counter,row_counter), &
               & constant)
       END DO
    END DO
    !$omp end do
    !$omp end parallel
  END SUBROUTINE ScaleDistributedSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the trace of the matrix.
  !! @param[in] this the matrix to compute the norm of.
  !! @return the trace value of the full distributed sparse matrix.
  FUNCTION Trace(this) RESULT(trace_value)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(in) :: this
    REAL(NTREAL) :: trace_value
    !! Local data
    TYPE(TripletList_t) :: triplet_list
    !! Counters/Temporary
    INTEGER :: counter
    TYPE(SparseMatrix_t) :: merged_local_data

    !! Merge all the local data
    CALL MergeLocalBlocks(this, merged_local_data)

    !! Compute The Local Contribution
    trace_value = 0
    CALL MatrixToTripletList(merged_local_data,triplet_list)
    DO counter = 1, triplet_list%CurrentSize
       IF (this%start_row + triplet_list%data(counter)%index_row .EQ. &
            & this%start_column + triplet_list%data(counter)%index_column) THEN
          trace_value = trace_value + triplet_list%data(counter)%point_value
       END IF
    END DO

    !! Sum Among Process Slice
    CALL MPI_Allreduce(MPI_IN_PLACE, trace_value, 1, MPINTREAL, &
         & MPI_SUM, within_slice_comm,grid_error)

    CALL DestructSparseMatrix(merged_local_data)
  END FUNCTION Trace
END MODULE DistributedSparseMatrixAlgebraModule
