  !! Communication Helpers
  TYPE(ReduceHelper_t), DIMENSION(:), ALLOCATABLE :: row_helper
  TYPE(ReduceHelper_t), DIMENSION(:), ALLOCATABLE :: column_helper
  TYPE(ReduceHelper_t), DIMENSION(:, :), ALLOCATABLE :: slice_helper
  !! For Iterating Over Local Blocks
  INTEGER :: II, II2, II2_range
  INTEGER :: JJ, JJ2, JJ2_range
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
  !! Temporary AB matrix for scaling.
  TYPE(Matrix_ps) :: matAB

  !! The threshold needs to be smaller if we are doing a sliced version
  !! because you might flush a value that would be kept in the summed version.
  IF (matA%process_grid%num_process_slices .GT. 1) THEN
     working_threshold = threshold / (matA%process_grid%num_process_slices*1000)
  ELSE
     working_threshold = threshold
  END IF

  !! Construct The Temporary Matrices
  CALL ConstructEmptyMatrix(matAB, matA)

  ALLOCATE(AdjacentABlocks(matAB%process_grid%number_of_blocks_rows, &
       & matAB%process_grid%number_of_blocks_columns / &
       & matAB%process_grid%num_process_slices))
  ALLOCATE(LocalRowContribution(matAB%process_grid%number_of_blocks_rows))
  ALLOCATE(GatheredRowContribution(matAB%process_grid%number_of_blocks_rows))
  ALLOCATE(GatheredRowContributionT(matAB%process_grid%number_of_blocks_rows))

  ALLOCATE(TransposedBBlocks(matAB%process_grid%number_of_blocks_rows / &
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
  DO II = 1, matAB%process_grid%number_of_blocks_rows
     ATasks(II) = LocalGatherA
  END DO
  ALLOCATE(BTasks(matAB%process_grid%number_of_blocks_columns))
  DO JJ = 1, matAB%process_grid%number_of_blocks_columns
     BTasks(JJ) = LocalGatherB
  END DO
  ALLOCATE(ABTasks(matAB%process_grid%number_of_blocks_rows, &
       & matAB%process_grid%number_of_blocks_columns))
  DO JJ = 1, matAB%process_grid%number_of_blocks_columns
     DO II = 1, matAB%process_grid%number_of_blocks_rows
        ABTasks(II,JJ) = AwaitingAB
     END DO
  END DO

  !! Setup A Tasks
  duplicate_start_column = matAB%process_grid%my_slice + 1
  duplicate_offset_column = matAB%process_grid%num_process_slices

  !! Setup B Tasks
  duplicate_start_row = matAB%process_grid%my_slice + 1
  duplicate_offset_row = matAB%process_grid%num_process_slices

  !! Run A Tasks
  ATasks_completed = 0
  BTasks_completed = 0
  ABTasks_completed = 0

  !$OMP PARALLEL
  !$OMP MASTER
  DO WHILE (ATasks_completed .LT. SIZE(ATasks) .OR. &
       & BTasks_completed .LT. SIZE(BTasks) .OR. &
       & ABTasks_completed .LT. SIZE(ABTasks))
     DO II = 1, matAB%process_grid%number_of_blocks_rows
        SELECT CASE (ATasks(II))
        CASE(LocalGatherA)
           ATasks(II) = TaskRunningA
           !$OMP TASK DEFAULT(SHARED), PRIVATE(JJ2, JJ2_range), FIRSTPRIVATE(II)
           !! First Align The Data We Are Working With
           JJ2_range = matAB%process_grid%number_of_blocks_columns / &
                     & matAB%process_grid%num_process_slices
           DO JJ2 = 1, JJ2_range
              CALL CopyMatrix(matA%LMAT(II, &
                   & duplicate_start_column + &
                   & duplicate_offset_column * (JJ2 - 1)),&
                   & AdjacentABlocks(II, JJ2))
           END DO
           !! Then Do A Local Gather and Cleanup
           CALL ComposeMatrixColumns(AdjacentABlocks(II, :), &
                & LocalRowContribution(II))
           DO JJ2 = 1, JJ2_range
              CALL DestructMatrix(AdjacentABlocks(II, JJ2))
           END DO
           ATasks(II) = SendSizeA
           !$OMP END TASK
        CASE(SendSizeA)
           !! Then Start A Global Gather
           CALL ReduceAndComposeMatrixSizes(LocalRowContribution(II), &
                & matAB%process_grid%blocked_row_comm(II), &
                & GatheredRowContribution(II), row_helper(II))
           ATasks(II) = ComposeA
        CASE(ComposeA)
           IF (TestReduceSizeRequest(row_helper(II))) THEN
              CALL ReduceAndComposeMatrixData(LocalRowContribution(II), &
                   & matAB%process_grid%blocked_row_comm(II), &
                   & GatheredRowContribution(II), row_helper(II))
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
           CALL ReduceAndComposeMatrixCleanup(LocalRowContribution(II), &
                & GatheredRowContribution(II), row_helper(II))
           CALL DestructMatrix(LocalRowContribution(II))
           CALL TransposeMatrix(GatheredRowContribution(II), &
                & GatheredRowContributionT(II))
           CALL DestructMatrix(GatheredRowContribution(II))
           ATasks(II) = CleanupA
           !$OMP END TASK
        CASE(CleanupA)
           ATasks(II) = FinishedA
           ATasks_completed = ATasks_completed + 1
        END SELECT
     END DO
     !! B Tasks
     DO JJ = 1 , matAB%process_grid%number_of_blocks_columns
        SELECT CASE (BTasks(JJ))
        CASE(LocalGatherB)
           BTasks(JJ) = TaskRunningB
           !$OMP TASK DEFAULT(SHARED), PRIVATE(II2, II2_range), FIRSTPRIVATE(JJ)
           !! First Transpose The Data We Are Working With
           II2_range = matAB%process_grid%number_of_blocks_rows / &
                     & matAB%process_grid%num_process_slices
           DO II2 = 1, II2_range
              CALL TransposeMatrix(matB%LMAT(duplicate_start_row + &
                   & duplicate_offset_row * (II2 - 1), JJ), &
                   & TransposedBBlocks(II2, JJ))
           END DO
           !! Then Do A Local Gather and Cleanup
           CALL ComposeMatrixColumns(TransposedBBlocks(:, JJ), &
                & LocalColumnContribution(JJ))
           DO II2 = 1, II2_range
              CALL DestructMatrix(TransposedBBlocks(II2, JJ))
           END DO
           BTasks(JJ) = SendSizeB
           !$OMP END TASK
        CASE(SendSizeB)
           !! Then A Global Gather
           CALL ReduceAndComposeMatrixSizes(LocalColumnContribution(JJ), &
                & matAB%process_grid%blocked_column_comm(JJ), &
                & GatheredColumnContribution(JJ), column_helper(JJ))
           BTasks(JJ) = LocalComposeB
        CASE(LocalComposeB)
           IF (TestReduceSizeRequest(column_helper(JJ))) THEN
              CALL ReduceAndComposeMatrixData(LocalColumnContribution(JJ),&
                   & matAB%process_grid%blocked_column_comm(JJ), &
                   & GatheredColumnContribution(JJ), column_helper(JJ))
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
           CALL ReduceAndComposeMatrixCleanup(LocalColumnContribution(JJ), &
                & GatheredColumnContribution(JJ), column_helper(JJ))
           CALL DestructMatrix(LocalColumnContribution(JJ))
           BTasks(JJ) = CleanupB
           !$OMP END TASK
        CASE(CleanupB)
           BTasks(JJ) = FinishedB
           BTasks_completed = BTasks_completed + 1
        END SELECT
     END DO
     !! AB Tasks
     DO II = 1 , matAB%process_grid%number_of_blocks_rows
        DO JJ = 1, matAB%process_grid%number_of_blocks_columns
           SELECT CASE(ABTasks(II, JJ))
           CASE (AwaitingAB)
              IF (ATasks(II) .EQ. FinishedA .AND. &
                   & BTasks(JJ) .EQ. FinishedB) THEN
                 ABTasks(II, JJ) = GemmAB
              END IF
           CASE (GemmAB)
              ABTasks(II, JJ) = TaskRunningAB
              !$OMP TASK DEFAULT(shared), FIRSTPRIVATE(II, JJ)
              CALL MatrixMultiply(GatheredRowContributionT(II), &
                   & GatheredColumnContribution(JJ), &
                   & SliceContribution(II, JJ), &
                   & IsATransposed_in = .TRUE., IsBTransposed_in = .TRUE., &
                   & alpha_in = alpha, threshold_in = working_threshold, &
                   & blocked_memory_pool_in = MPGRID(II, JJ))
              !! We can exit early if there is only one process slice
              IF (matAB%process_grid%num_process_slices .EQ. 1) THEN
                 ABTasks(II,JJ) = CleanupAB
                 CALL CopyMatrix(SliceContribution(II, JJ), matAB%LMAT(II, JJ))
                 CALL DestructMatrix(SliceContribution(II, JJ))
              ELSE
                 ABTasks(II, JJ) = SendSizeAB
              END IF
              !$OMP END TASK
           CASE(SendSizeAB)
              CALL ReduceAndSumMatrixSizes(SliceContribution(II, JJ),&
                   & matAB%process_grid%blocked_between_slice_comm(II, JJ), &
                   & matAB%LMAT(II, JJ), slice_helper(II, JJ))
              ABTasks(II, JJ) = GatherAndSumAB
           CASE (GatherAndSumAB)
              IF (TestReduceSizeRequest(slice_helper(II, JJ))) THEN
                 CALL ReduceAndSumMatrixData(SliceContribution(II, JJ), &
                      & matAB%process_grid%blocked_between_slice_comm(II, JJ), &
                      & matAB%LMAT(II, JJ), slice_helper(II, JJ))
                 ABTasks(II, JJ) = WaitInnerAB
              END IF
           CASE (WaitInnerAB)
              IF (TestReduceInnerRequest(slice_helper(II, JJ))) THEN
                 ABTasks(II, JJ) = WaitDataAB
              END IF
           CASE (WaitDataAB)
              IF (TestReduceDataRequest(slice_helper(II, JJ))) THEN
                 ABTasks(II, JJ) = LocalSumAB
              END IF
           CASE(LocalSumAB)
              ABTasks(II, JJ) = TaskRunningAB
              !$OMP TASK DEFAULT(SHARED), FIRSTPRIVATE(II, JJ)
              CALL ReduceAndSumMatrixCleanup(SliceContribution(II, JJ), &
                   & matAB%LMAT(II, JJ), threshold, slice_helper(II, JJ))
              CALL DestructMatrix(SliceContribution(II, JJ))
              ABTasks(II, JJ) = CleanupAB
              !$OMP END TASK
           CASE(CleanupAB)
              ABTasks(II, JJ) = FinishedAB
              ABTasks_completed = ABTasks_completed + 1
           END SELECT
        END DO
     END DO
     !! Prevent deadlock in the case where the number of tasks is capped.
     IF (matA%process_grid%omp_max_threads .EQ. 1) THEN
        !$OMP taskwait
     END IF
  END DO
  !$OMP END MASTER
  !$OMP END PARALLEL

  !! Cleanup
  DEALLOCATE(row_helper)
  DEALLOCATE(column_helper)
  DEALLOCATE(slice_helper)
  DEALLOCATE(ATasks)
  DEALLOCATE(BTasks)
  DEALLOCATE(ABTasks)

  !! Deallocate Buffers From A
  DO II = 1, matAB%process_grid%number_of_blocks_rows
     DO JJ2 = 1, matAB%process_grid%number_of_blocks_columns / &
          & matAB%process_grid%num_process_slices
        CALL DestructMatrix(AdjacentABlocks(II, JJ2))
     END DO
     CALL DestructMatrix(LocalRowContribution(II))
     CALL DestructMatrix(GatheredRowContribution(II))
  END DO
  DEALLOCATE(AdjacentABlocks)
  DEALLOCATE(LocalRowContribution)
  DEALLOCATE(GatheredRowContribution)
  !! Deallocate Buffers From B
  DO JJ = 1, matAB%process_grid%number_of_blocks_columns
     DO II2 = 1, matAB%process_grid%number_of_blocks_rows / &
          & matAB%process_grid%num_process_slices
        CALL DestructMatrix(TransposedBBlocks(II2, JJ))
     END DO
     CALL DestructMatrix(LocalColumnContribution(JJ))
  END DO
  DEALLOCATE(TransposedBBlocks)
  DEALLOCATE(LocalColumnContribution)
  !! Deallocate Buffers From Multiplying The Block
  DO II = 1, matAB%process_grid%number_of_blocks_rows
     CALL DestructMatrix(GatheredRowContributionT(II))
  END DO
  DO JJ = 1, matAB%process_grid%number_of_blocks_columns
     CALL DestructMatrix(GatheredColumnContribution(JJ))
  END DO
  DEALLOCATE(GatheredRowContributionT)
  DEALLOCATE(GatheredColumnContribution)
  !! Deallocate Buffers From Sum
  DO JJ = 1, matAB%process_grid%number_of_blocks_columns
     DO II = 1, matAB%process_grid%number_of_blocks_rows
        CALL DestructMatrix(SliceContribution(II, JJ))
     END DO
  END DO
  DEALLOCATE(SliceContribution)

  !! Copy to output matrix.
  IF (ABS(beta) .LT. TINY(beta)) THEN
     CALL CopyMatrix(matAB, matC)
  ELSE
     CALL ScaleMatrix(MatC, beta)
     CALL IncrementMatrix(MatAB, MatC)
  END IF
  CALL DestructMatrix(matAB)
