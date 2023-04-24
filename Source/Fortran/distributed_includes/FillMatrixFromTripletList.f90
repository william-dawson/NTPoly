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

  IF (prepartitioned) THEN
     !! Shift and sort the local entries.
     shifted = triplet_list
     CALL ShiftTripletList(shifted, 1 - this%start_row, 1 - this%start_column)
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

     !! And reduce over the Z dimension. 
     IF (.NOT. preduplicated .AND. &
          & .NOT. this%process_grid%num_process_slices .EQ. 1) THEN
        CALL ReduceAndSumMatrix(local_matrix, gathered_matrix, threshold, &
             & this%process_grid%between_slice_comm)
        CALL SplitMatrixToLocalBlocks(this, gathered_matrix)
     ELSE
        CALL SplitMatrixToLocalBlocks(this, local_matrix)
     END IF
  END IF

  CALL DestructMatrix(local_matrix)
  CALL DestructTripletList(sorted_triplet_list)
