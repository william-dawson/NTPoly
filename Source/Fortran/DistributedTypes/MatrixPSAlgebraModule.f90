!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Performing Distributed Sparse Matrix Algebra Operations.
MODULE MatrixPSAlgebraModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL
  USE MatrixMemoryPoolPModule, ONLY : MatrixMemoryPool_p, &
       & CheckMemoryPoolValidity, ConstructMatrixMemoryPool
  USE MatrixPSModule
  USE GemmTasksModule
  USE MatrixReduceModule, ONLY : ReduceHelper_t, ReduceSizes, &
       & ReduceAndComposeData, ReduceAndComposeCleanup, ReduceAndSumData, &
       & ReduceAndSumCleanup, TestReduceSizeRequest, TestReduceOuterRequest, &
       & TestReduceInnerRequest, TestReduceDataRequest
  USE MatrixSAlgebraModule
  USE MatrixSModule, ONLY : Matrix_lsr, DestructMatrix,  CopyMatrix, &
       & TransposeMatrix, ComposeMatrixColumns, MatrixToTripletList
  USE TimerModule, ONLY : StartTimer, StopTimer
  USE TripletListModule, ONLY : TripletList_r, DestructTripletList
  USE ISO_C_BINDING
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ComputeMatrixSigma
  PUBLIC :: MatrixMultiply
  PUBLIC :: MatrixGrandSum
  PUBLIC :: PairwiseMultiplyMatrix
  PUBLIC :: MatrixNorm
  PUBLIC :: DotMatrix
  PUBLIC :: IncrementMatrix
  PUBLIC :: ScaleMatrix
  PUBLIC :: MatrixTrace
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE ComputeMatrixSigma
    MODULE PROCEDURE ComputeMatrixSigma_ps
  END INTERFACE
  INTERFACE MatrixMultiply
    MODULE PROCEDURE MatrixMultiply_ps
  END INTERFACE
  INTERFACE MatrixGrandSum
    MODULE PROCEDURE MatrixGrandSum_ps
  END INTERFACE
  INTERFACE PairwiseMultiplyMatrix
    MODULE PROCEDURE PairwiseMultiplyMatrix_ps
  END INTERFACE
  INTERFACE MatrixNorm
    MODULE PROCEDURE MatrixNorm_ps
  END INTERFACE
  INTERFACE DotMatrix
    MODULE PROCEDURE DotMatrix_ps
  END INTERFACE
  INTERFACE IncrementMatrix
    MODULE PROCEDURE IncrementMatrix_ps
  END INTERFACE
  INTERFACE ScaleMatrix
    MODULE PROCEDURE ScaleMatrix_ps
  END INTERFACE
  INTERFACE MatrixTrace
    MODULE PROCEDURE MatrixTrace_ps
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute sigma for the inversion method.
  !! See \cite ozaki2001efficient for details.
  !! @param[in] this the matrix to compute the sigma value of.
  !! @param[out] sigma_value sigma.
  SUBROUTINE ComputeMatrixSigma_ps(this, sigma_value)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: this
    REAL(NTREAL), INTENT(OUT) :: sigma_value
    !! Local Data
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: column_sigma_contribution
    !! Counters/Temporary
    INTEGER :: inner_counter, outer_counter
    TYPE(Matrix_lsr) :: merged_local_data
    INTEGER :: ierr

    !! Merge all the local data
    CALL MergeMatrixLocalBlocks(this, merged_local_data)

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
         & this%process_grid%column_comm, ierr)
    CALL MPI_Allreduce(MAXVAL(column_sigma_contribution),sigma_value,1, &
         & MPINTREAL,MPI_MAX, this%process_grid%row_comm, ierr)
    sigma_value = 1.0d+0/(sigma_value**2)

    DEALLOCATE(column_sigma_contribution)
    CALL DestructMatrix(merged_local_data)
  END SUBROUTINE ComputeMatrixSigma_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Multiply two matrices together, and add to the third.
  !! C := alpha*matA*matB+ beta*matC
  !! @param[in] matA Matrix A
  !! @param[in] matB Matrix B
  !! @param[out] matC = alpha*matA*matB + beta*matC
  !! @param[in] alpha_in scales the multiplication
  !! @param[in] beta_in scales matrix we sum on to
  !! @param[in] threshold_in for flushing values to zero (Optional, default=0).
  !! @param[inout] memory_pool_in a memory pool for the calculation (Optional).
  SUBROUTINE MatrixMultiply_ps(matA, matB ,matC, alpha_in, beta_in, threshold_in,&
       & memory_pool_in)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN)    :: matA
    TYPE(Matrix_ps), INTENT(IN)    :: matB
    TYPE(Matrix_ps), INTENT(INOUT) :: matC
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: beta_in
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    TYPE(MatrixMemoryPool_p), OPTIONAL, INTENT(inout) :: &
         & memory_pool_in
    !! Local Versions of Optional Parameter
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: beta
    REAL(NTREAL) :: threshold
    TYPE(Matrix_ps) :: matAB
    !! Temporary Matrices
    TYPE(Matrix_lsr), DIMENSION(:,:), ALLOCATABLE :: AdjacentABlocks
    TYPE(Matrix_lsr), DIMENSION(:), ALLOCATABLE :: LocalRowContribution
    TYPE(Matrix_lsr), DIMENSION(:), ALLOCATABLE :: GatheredRowContribution
    TYPE(Matrix_lsr), DIMENSION(:), ALLOCATABLE :: GatheredRowContributionT
    TYPE(Matrix_lsr), DIMENSION(:,:), ALLOCATABLE :: TransposedBBlocks
    TYPE(Matrix_lsr), DIMENSION(:), ALLOCATABLE :: LocalColumnContribution
    TYPE(Matrix_lsr), DIMENSION(:), ALLOCATABLE :: GatheredColumnContribution
    TYPE(Matrix_lsr), DIMENSION(:,:), ALLOCATABLE :: SliceContribution
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
    IF (matA%process_grid%num_process_slices .GT. 1) THEN
       working_threshold = threshold/(matA%process_grid%num_process_slices*1000)
    ELSE
       working_threshold = threshold
    END IF

    !! Construct The Temporary Matrices
    CALL ConstructEmptyMatrix(matAB, &
         & matA%actual_matrix_dimension, matA%process_grid)

    ALLOCATE(AdjacentABlocks(matAB%process_grid%number_of_blocks_rows, &
         & matAB%process_grid%number_of_blocks_columns/&
         & matAB%process_grid%num_process_slices))
    ALLOCATE(LocalRowContribution(matAB%process_grid%number_of_blocks_rows))
    ALLOCATE(GatheredRowContribution(matAB%process_grid%number_of_blocks_rows))
    ALLOCATE(GatheredRowContributionT(matAB%process_grid%number_of_blocks_rows))

    ALLOCATE(TransposedBBlocks(matAB%process_grid%number_of_blocks_rows/&
         & matAB%process_grid%num_process_slices, &
         & matAB%process_grid%number_of_blocks_columns))
    ALLOCATE(LocalColumnContribution(&
         & matAB%process_grid%number_of_blocks_columns))
    ALLOCATE(GatheredColumnContribution(&
         & matAB%process_grid%number_of_blocks_columns))
    ALLOCATE(SliceContribution(matAB%process_grid%number_of_blocks_rows, &
         & matAB%process_grid%number_of_blocks_columns))

    !! Helpers
    ALLOCATE(row_helper(matAB%process_grid%number_of_blocks_rows))
    ALLOCATE(column_helper(matAB%process_grid%number_of_blocks_columns))
    ALLOCATE(slice_helper(matAB%process_grid%number_of_blocks_rows, &
         & matAB%process_grid%number_of_blocks_columns))

    !! Construct the task queues
    ALLOCATE(ATasks(matAB%process_grid%number_of_blocks_rows))
    DO II=1,matAB%process_grid%number_of_blocks_rows
       ATasks(II) = LocalGatherA
    END DO
    ALLOCATE(BTasks(matAB%process_grid%number_of_blocks_columns))
    DO JJ=1,matAB%process_grid%number_of_blocks_columns
       BTasks(JJ) = LocalGatherB
    END DO
    ALLOCATE(ABTasks(matAB%process_grid%number_of_blocks_rows, &
         & matAB%process_grid%number_of_blocks_columns))
    DO JJ=1,matAB%process_grid%number_of_blocks_columns
       DO II=1,matAB%process_grid%number_of_blocks_rows
          ABTasks(II,JJ) = AwaitingAB
       END DO
    END DO

    !! Setup A Tasks
    duplicate_start_column = matAB%process_grid%my_slice+1
    duplicate_offset_column = matAB%process_grid%num_process_slices

    !! Setup B Tasks
    duplicate_start_row = matAB%process_grid%my_slice+1
    duplicate_offset_row = matAB%process_grid%num_process_slices

    !! Setup AB Tasks
    IF (PRESENT(memory_pool_in)) THEN
       IF (.NOT. CheckMemoryPoolValidity(memory_pool_in)) THEN
          CALL ConstructMatrixMemoryPool(memory_pool_in, matAB)
       END IF
    END IF

    !! Run A Tasks
    ATasks_completed = 0
    BTasks_completed = 0
    ABTasks_completed = 0
    !$OMP PARALLEL
    !$OMP MASTER
    DO WHILE (ATasks_completed .LT. SIZE(ATasks) .OR. &
         & BTasks_completed .LT. SIZE(BTasks) .OR. &
         & ABTasks_completed .LT. SIZE(ABTasks))
       DO II=1, matAB%process_grid%number_of_blocks_rows
          SELECT CASE (ATasks(II))
          CASE(LocalGatherA)
             ATasks(II) = TaskRunningA
             !$OMP TASK DEFAULT(SHARED), PRIVATE(JJ2), FIRSTPRIVATE(II)
             !! First Align The Data We're Working With
             DO JJ2=1, &
                  & matAB%process_grid%number_of_blocks_columns/ &
                  & matAB%process_grid%num_process_slices
                CALL CopyMatrix(matA%local_data(II, &
                     & duplicate_start_column+duplicate_offset_column*(JJ2-1)),&
                     & AdjacentABlocks(II,JJ2))
             END DO
             !! Then Do A Local Gather
             CALL ComposeMatrixColumns(AdjacentABlocks(II,:), &
                  & LocalRowContribution(II))
             ATasks(II) = SendSizeA
             !$OMP END TASK
          CASE(SendSizeA)
             !! Then Start A Global Gather
             CALL ReduceSizes(LocalRowContribution(II), &
                  & matAB%process_grid%blocked_row_comm(II), &
                  & row_helper(II))
             ATasks(II) = ComposeA
          CASE(ComposeA)
             IF (TestReduceSizeRequest(row_helper(II))) THEN
                CALL ReduceAndComposeData(LocalRowContribution(II), &
                     & matAB%process_grid%blocked_row_comm(II), &
                     & GatheredRowContribution(II), row_helper(II))
                ATasks(II) = WaitOuterA
             END IF
          CASE(WaitOuterA)
             IF (TestReduceOuterRequest(row_helper(II))) THEN
                ATasks(II) = WaitInnerA
             END IF
          CASE(WaitInnerA)
             IF (TestReduceInnerRequest(row_helper(II))) THEN
                ATasks(II) = WaitDataA
             END IF
          CASE(WaitDataA)
             IF (TestReduceDataRequest(row_helper(II))) THEN
                ATasks(II) = AdjustIndicesA
             END IF
          CASE(AdjustIndicesA)
             ATasks(II) = TaskRunningA
             !$OMP TASK DEFAULT(SHARED), FIRSTPRIVATE(II)
             CALL ReduceAndComposeCleanup(LocalRowContribution(II), &
                  & GatheredRowContribution(II), row_helper(II))
             CALL TransposeMatrix(GatheredRowContribution(II), &
                  & GatheredRowContributionT(II))
             ATasks(II) = CleanupA
             !$OMP END TASK
          CASE(CleanupA)
             ATasks(II) = FinishedA
             ATasks_completed = ATasks_completed + 1
          END SELECT
       END DO
       !! B Tasks
       DO JJ=1,matAB%process_grid%number_of_blocks_columns
          SELECT CASE (BTasks(JJ))
          CASE(LocalGatherB)
             BTasks(JJ) = TaskRunningB
             !$OMP TASK DEFAULT(SHARED), PRIVATE(II2), FIRSTPRIVATE(JJ)
             !! First Transpose The Data We're Working With
             DO II2=1, matAB%process_grid%number_of_blocks_rows/&
                  & matAB%process_grid%num_process_slices
                CALL TransposeMatrix(matB%local_data(duplicate_start_row+&
                     & duplicate_offset_row*(II2-1),JJ), &
                     & TransposedBBlocks(II2,JJ))
             END DO
             !! Then Do A Local Gather
             CALL ComposeMatrixColumns(TransposedBBlocks(:,JJ), &
                  & LocalColumnContribution(JJ))
             BTasks(JJ) = SendSizeB
             !$OMP END TASK
          CASE(SendSizeB)
             !! Then A Global Gather
             CALL ReduceSizes(LocalColumnContribution(JJ), &
                  & matAB%process_grid%blocked_column_comm(JJ), &
                  & column_helper(JJ))
             BTasks(JJ) = LocalComposeB
          CASE(LocalComposeB)
             IF (TestReduceSizeRequest(column_helper(JJ))) THEN
                CALL ReduceAndComposeData(LocalColumnContribution(JJ),&
                     & matAB%process_grid%blocked_column_comm(JJ), &
                     & GatheredColumnContribution(JJ), column_helper(JJ))
                BTasks(JJ) = WaitOuterB
             END IF
          CASE(WaitOuterB)
             IF (TestReduceOuterRequest(column_helper(JJ))) THEN
                BTasks(JJ) = WaitInnerB
             END IF
          CASE(WaitInnerB)
             IF (TestReduceInnerRequest(column_helper(JJ))) THEN
                BTasks(JJ) = WaitDataB
             END IF
          CASE(WaitDataB)
             IF (TestReduceDataRequest(column_helper(JJ))) THEN
                BTasks(JJ) = AdjustIndicesB
             END IF
          CASE(AdjustIndicesB)
             BTasks(JJ) = TaskRunningB
             !$OMP TASK DEFAULT(SHARED), FIRSTPRIVATE(JJ)
             CALL ReduceAndComposeCleanup(LocalColumnContribution(JJ), &
                  & GatheredColumnContribution(JJ), column_helper(JJ))
             BTasks(JJ) = CleanupB
             !$OMP END TASK
          CASE(CleanupB)
             BTasks(JJ) = FinishedB
             BTasks_completed = BTasks_completed + 1
          END SELECT
       END DO
       !! AB Tasks
       DO II=1,matAB%process_grid%number_of_blocks_rows
          DO JJ=1,matAB%process_grid%number_of_blocks_columns
             SELECT CASE(ABTasks(II,JJ))
             CASE (AwaitingAB)
                IF (ATasks(II) .EQ. FinishedA .AND. &
                     & BTasks(JJ) .EQ. FinishedB) THEN
                   ABTasks(II,JJ) = GemmAB
                END IF
             CASE (GemmAB)
                ABTasks(II,JJ) = TaskRunningAB
                !$OMP TASK DEFAULT(shared), FIRSTPRIVATE(II,JJ)
                IF (PRESENT(memory_pool_in)) THEN
                   CALL MatrixMultiply(GatheredRowContributionT(II), &
                        & GatheredColumnContribution(JJ), &
                        & SliceContribution(II,JJ), &
                        & IsATransposed_in=.TRUE., IsBTransposed_in=.TRUE., &
                        & alpha_in=alpha, threshold_in=working_threshold, &
                        & blocked_memory_pool_in=memory_pool_in%grid(II,JJ))
                ELSE
                   CALL MatrixMultiply(GatheredRowContributionT(II), &
                        & GatheredColumnContribution(JJ), &
                        & SliceContribution(II,JJ), &
                        & IsATransposed_in=.TRUE., IsBTransposed_in=.TRUE., &
                        & alpha_in=alpha, threshold_in=working_threshold)
                END IF
                ABTasks(II,JJ) = SendSizeAB
                !$OMP END TASK
             CASE(SendSizeAB)
                CALL ReduceSizes(SliceContribution(II,JJ),&
                     & matAB%process_grid%blocked_between_slice_comm(II,JJ), &
                     & slice_helper(II,JJ))
                ABTasks(II,JJ) = GatherAndSumAB
             CASE (GatherAndSumAB)
                IF (TestReduceSizeRequest(&
                     & slice_helper(II,JJ))) THEN
                   CALL ReduceAndSumData(SliceContribution(II,JJ), &
                        & matAB%process_grid%blocked_between_slice_comm(II,JJ),&
                        & slice_helper(II,JJ))
                   ABTasks(II,JJ) = WaitOuterAB
                END IF
             CASE (WaitOuterAB)
                IF (TestReduceOuterRequest(&
                     & slice_helper(II,JJ))) THEN
                   ABTasks(II,JJ) = WaitInnerAB
                END IF
             CASE (WaitInnerAB)
                IF (TestReduceInnerRequest(&
                     & slice_helper(II,JJ))) THEN
                   ABTasks(II,JJ) = WaitDataAB
                END IF
             CASE (WaitDataAB)
                IF (TestReduceDataRequest(&
                     & slice_helper(II,JJ))) THEN
                   ABTasks(II,JJ) = LocalSumAB
                END IF
             CASE(LocalSumAB)
                ABTasks(II,JJ) = TaskRunningAB
                !$OMP TASK DEFAULT(SHARED), FIRSTPRIVATE(II,JJ)
                CALL ReduceAndSumCleanup(SliceContribution(II,JJ), &
                     & matAB%local_data(II,JJ), threshold, slice_helper(II,JJ))
                ABTasks(II,JJ) = CleanupAB
                !$OMP END TASK
             CASE(CleanupAB)
                ABTasks(II,JJ) = FinishedAB
                ABTasks_completed = ABTasks_completed + 1
             END SELECT
          END DO
       END DO
    END DO
    !$OMP END MASTER
    !$OMP END PARALLEL

    !! Copy to output matrix.
    IF (beta .EQ. 0.0) THEN
       CALL CopyMatrix(matAB,matC)
    ELSE
       CALL ScaleMatrix(MatC,beta)
       CALL IncrementMatrix(MatAB,MatC)
    END IF

    !! Cleanup
    CALL DestructMatrix(matAB)
    DEALLOCATE(row_helper)
    DEALLOCATE(column_helper)
    DEALLOCATE(slice_helper)

    !! Deallocate Buffers From A
    DO II=1,matAB%process_grid%number_of_blocks_rows
       DO JJ2=1,matAB%process_grid%number_of_blocks_columns/&
            & matAB%process_grid%num_process_slices
          CALL DestructMatrix(AdjacentABlocks(II,JJ2))
       END DO
       CALL DestructMatrix(LocalRowContribution(II))
       CALL DestructMatrix(GatheredRowContribution(II))
    END DO
    DEALLOCATE(AdjacentABlocks)
    DEALLOCATE(LocalRowContribution)
    DEALLOCATE(GatheredRowContribution)
    !! Deallocate Buffers From B
    DO JJ=1,matAB%process_grid%number_of_blocks_columns
       DO II2=1,matAB%process_grid%number_of_blocks_rows/&
            & matAB%process_grid%num_process_slices
          CALL DestructMatrix(TransposedBBlocks(II2,JJ))
       END DO
       CALL DestructMatrix(LocalColumnContribution(JJ))
    END DO
    DEALLOCATE(TransposedBBlocks)
    DEALLOCATE(LocalColumnContribution)
    !! Deallocate Buffers From Multiplying The Block
    DO II=1,matAB%process_grid%number_of_blocks_rows
       CALL DestructMatrix(GatheredRowContributionT(II))
    END DO
    DO JJ=1,matAB%process_grid%number_of_blocks_columns
       CALL DestructMatrix(GatheredColumnContribution(JJ))
    END DO
    DEALLOCATE(GatheredRowContributionT)
    DEALLOCATE(GatheredColumnContribution)
    !! Deallocate Buffers From Sum
    DO JJ=1,matAB%process_grid%number_of_blocks_columns
       DO II=1,matAB%process_grid%number_of_blocks_rows
          CALL DestructMatrix(SliceContribution(II,JJ))
       END DO
    END DO
    DEALLOCATE(SliceContribution)

    CALL StopTimer("GEMM")
  END SUBROUTINE MatrixMultiply_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sum up the elements in a matrix into a single value.
  !! @param[in] this matrix to compute.
  !! @result sum the sum of all elements.
  FUNCTION MatrixGrandSum_ps(this) RESULT(sum)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN)  :: this
    REAL(NTREAL) :: sum
    !! Local Data
    INTEGER :: II, JJ
    REAL(NTREAL) :: temp
    INTEGER :: ierr

    sum = 0
    DO JJ = 1, this%process_grid%number_of_blocks_columns
       DO II = 1, this%process_grid%number_of_blocks_rows
          temp = MatrixGrandSum(this%local_data(II,JJ))
          sum = sum + temp
       END DO
    END DO

    !! Sum Among Process Slice
    CALL MPI_Allreduce(MPI_IN_PLACE, sum, 1, MPINTREAL, &
         & MPI_SUM, this%process_grid%within_slice_comm, ierr)
  END FUNCTION MatrixGrandSum_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Elementwise multiplication. C_ij = A_ij * B_ij.
  !! Also known as a Hadamard product.
  !! @param[in] matA Matrix A.
  !! @param[in] matB Matrix B.
  !! @param[inout] matC = MatA mult MatB.
  SUBROUTINE PairwiseMultiplyMatrix_ps(matA, matB, matC)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN)  :: matA
    TYPE(Matrix_ps), INTENT(IN)  :: matB
    TYPE(Matrix_ps), INTENT(INOUT)  :: matC
    !! Local Data
    INTEGER :: II, JJ

    CALL ConstructEmptyMatrix(matC, &
         & matA%actual_matrix_dimension, matA%process_grid)

    !$omp parallel
    !$omp do collapse(2)
    DO JJ = 1, matA%process_grid%number_of_blocks_columns
       DO II = 1, matA%process_grid%number_of_blocks_rows
          CALL PairwiseMultiplyMatrix(matA%local_data(II,JJ), &
               & matB%local_data(II,JJ), matC%local_data(II,JJ))
       END DO
    END DO
    !$omp end do
    !$omp end parallel
  END SUBROUTINE PairwiseMultiplyMatrix_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the norm of a distributed sparse matrix along the rows.
  !! @param[in] this the matrix to compute the norm of.
  !! @return the norm value of the full distributed sparse matrix.
  FUNCTION MatrixNorm_ps(this) RESULT(norm_value)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: this
    REAL(NTREAL) :: norm_value
    !! Local Data
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: local_norm
    TYPE(Matrix_lsr) :: merged_local_data
    INTEGER :: ierr

    !! Merge all the local data
    CALL MergeMatrixLocalBlocks(this, merged_local_data)
    ALLOCATE(local_norm(merged_local_data%columns))

    !! Sum Along Columns
    CALL MatrixColumnNorm(merged_local_data,local_norm)
    CALL MPI_Allreduce(MPI_IN_PLACE,local_norm,SIZE(local_norm), &
         & MPINTREAL, MPI_SUM, this%process_grid%column_comm, ierr)

    !! Find Max Value Amonst Columns
    norm_value = MAXVAL(local_norm)
    CALL MPI_Allreduce(MPI_IN_PLACE,norm_value,1,MPINTREAL,MPI_MAX, &
         & this%process_grid%row_comm, ierr)

    CALL DestructMatrix(merged_local_data)
    DEALLOCATE(local_norm)
  END FUNCTION MatrixNorm_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> product = dot(Matrix A,Matrix B)
  !! Note that a dot product is the sum of elementwise multiplication, not
  !! traditional matrix multiplication.
  !! @param[in] matA Matrix A.
  !! @param[in,out] matB Matrix B.
  !! @result product the dot product.
  FUNCTION DotMatrix_ps(matA, matB) RESULT(product)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN)  :: matA
    TYPE(Matrix_ps), INTENT(IN)  :: matB
    REAL(NTREAL) :: product
    !! Local Data
    TYPE(Matrix_ps)  :: matC

    CALL PairwiseMultiplyMatrix(matA,matB,matC)
    product = MatrixGrandSum(matC)
    CALL DestructMatrix(matC)
  END FUNCTION DotMatrix_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Matrix B = alpha*Matrix A + Matrix B (AXPY)
  !! This will utilize the sparse vector increment routine.
  !! @param[in] matA Matrix A.
  !! @param[inout] matB Matrix B.
  !! @param[in] alpha_in multiplier. (Optional, default= 1.0)
  !! @param[in] threshold_in for flushing values to zero. (Optional, default=0)
  SUBROUTINE IncrementMatrix_ps(matA, matB, alpha_in,threshold_in)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN)  :: matA
    TYPE(Matrix_ps), INTENT(INOUT)  :: matB
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    !! Local Data
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

    !$omp parallel
    !$omp do collapse(2)
    DO JJ = 1, matA%process_grid%number_of_blocks_columns
       DO II = 1, matA%process_grid%number_of_blocks_rows
          CALL IncrementMatrix(matA%local_data(II,JJ), &
               & matB%local_data(II,JJ), alpha, threshold)
       END DO
    END DO
    !$omp end do
    !$omp end parallel

  END SUBROUTINE IncrementMatrix_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a distributed sparse matrix by a constant.
  !! @param[inout] this Matrix to scale.
  !! @param[in] constant scale factor.
  SUBROUTINE ScaleMatrix_ps(this,constant)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    REAL(NTREAL), INTENT(IN) :: constant
    INTEGER :: II, JJ

    !$omp parallel
    !$omp do collapse(2)
    DO JJ = 1, this%process_grid%number_of_blocks_columns
       DO II = 1, this%process_grid%number_of_blocks_rows
          CALL ScaleMatrix(this%local_data(II,JJ),constant)
       END DO
    END DO
    !$omp end do
    !$omp end parallel
  END SUBROUTINE ScaleMatrix_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the trace of the matrix.
  !! @param[in] this the matrix to compute the trace of.
  !! @return the trace value of the full distributed sparse matrix.
  FUNCTION MatrixTrace_ps(this) RESULT(trace_value)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: this
    REAL(NTREAL) :: trace_value
    !! Local data
    TYPE(TripletList_r) :: triplet_list
    !! Counters/Temporary
    INTEGER :: counter
    TYPE(Matrix_lsr) :: merged_local_data
    INTEGER :: ierr

    !! Merge all the local data
    CALL MergeMatrixLocalBlocks(this, merged_local_data)

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
         & MPI_SUM, this%process_grid%within_slice_comm, ierr)

    CALL DestructMatrix(merged_local_data)
  END FUNCTION MatrixTrace_ps
END MODULE MatrixPSAlgebraModule
