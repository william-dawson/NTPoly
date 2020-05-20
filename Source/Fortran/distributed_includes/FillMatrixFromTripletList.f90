  !! Optional Parameteres
  IF (.NOT. PRESENT(preduplicated_in)) THEN
     preduplicated = .FALSE.
  ELSE
     preduplicated = preduplicated_in
  END IF

  IF (.NOT. PRESENT(prepartitioned_in)) THEN
     prepartitioned = .FALSE.
  ELSE
     prepartitioned = prepartitioned_in
  END IF

  CALL StartTimer("FillFromTriplet")

  IF (prepartitioned) THEN
     !! Shift and sort the local entries.
     CALL ConstructTripletList(shifted, triplet_list%CurrentSize)
     shifted%data(:triplet_list%CurrentSize) = triplet_list%data(:triplet_list%CurrentSize)
     DO II = 1, shifted%CurrentSize
        local_column = shifted%data(II)%index_column - this%start_column + 1
        local_row = shifted%data(II)%index_row - this%start_row + 1
        shifted%data(II)%index_column = local_column
        shifted%data(II)%index_row = local_row
     END DO
     CALL SortTripletList(shifted, this%local_columns, &
          & this%local_rows, sorted_triplet_list)
     !! Build
     CALL ConstructMatrixFromTripletList(local_matrix, sorted_triplet_list, &
          & this%local_rows, this%local_columns)
     CALL SplitMatrixToLocalBlocks(this, local_matrix)
  ELSE
     !! First we redistribute the triplet list to get all the local data
     !! on the correct process.
     CALL ConstructDefaultPermutation(basic_permutation, &
          & this%logical_matrix_dimension)
     CALL RedistributeData(this,basic_permutation%index_lookup, &
          & basic_permutation%reverse_index_lookup, triplet_list, &
          & sorted_triplet_list)

     !! Now we can just construct a local matrix.
     CALL ConstructMatrixFromTripletList(local_matrix, sorted_triplet_list, &
          & this%local_rows, this%local_columns)
     !! And reduce over the Z dimension. This can be accomplished by
     !! summing up.
     IF (.NOT. preduplicated .AND. &
          & .NOT. this%process_grid%num_process_slices .EQ. 1) THEN
        CALL ReduceAndSumMatrixSizes(local_matrix, &
             & this%process_grid%between_slice_comm, gathered_matrix, &
             & gather_helper)
        DO WHILE(.NOT. TestReduceSizeRequest(gather_helper))
        END DO
        CALL ReduceAndSumMatrixData(local_matrix, gathered_matrix, &
             & this%process_grid%between_slice_comm, gather_helper)
        DO WHILE(.NOT. TestReduceInnerRequest(gather_helper))
        END DO
        DO WHILE(.NOT. TestReduceDataRequest(gather_helper))
        END DO
        CALL ReduceAndSumMatrixCleanup(local_matrix, gathered_matrix, &
             & threshold, gather_helper)
        CALL SplitMatrixToLocalBlocks(this, gathered_matrix)
     ELSE
        CALL SplitMatrixToLocalBlocks(this, local_matrix)
     END IF
     CALL StopTimer("FillFromTriplet")
  END IF
  CALL DestructMatrix(local_matrix)
  CALL DestructTripletList(sorted_triplet_list)
