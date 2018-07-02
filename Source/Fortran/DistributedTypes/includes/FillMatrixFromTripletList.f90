  !! Local Data
  TYPE(Permutation_t) :: basic_permutation
  TYPE(ReduceHelper_t) :: gather_helper
  REAL(NTREAL), PARAMETER :: threshold = 0.0
  LOGICAL :: preduplicated

  IF (.NOT. PRESENT(preduplicated_in)) THEN
     preduplicated = .FALSE.
  ELSE
     preduplicated = preduplicated_in
  END IF

  CALL StartTimer("FillFromTriplet")
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
  IF (.NOT. PRESENT(preduplicated_in) .OR. .NOT. preduplicated_in) THEN
     CALL ReduceMatrixSizes(local_matrix, this%process_grid%between_slice_comm,&
          & gather_helper)
     DO WHILE(.NOT. TestReduceSizeRequest(gather_helper))
     END DO
     CALL ReduceAndSumMatrixData(local_matrix, gathered_matrix, &
          & this%process_grid%between_slice_comm, gather_helper)
     DO WHILE(.NOT. TestReduceOuterRequest(gather_helper))
     END DO
     DO WHILE(.NOT. TestReduceInnerRequest(gather_helper))
     END DO
     DO WHILE(.NOT. TestReduceDataRequest(gather_helper))
     END DO
     CALL ReduceAndSumMatrixCleanup(local_matrix, gathered_matrix, threshold, &
          & gather_helper)
     CALL SplitMatrixToLocalBlocks(this, gathered_matrix)
  ELSE
     CALL SplitMatrixToLocalBlocks(this, local_matrix)
  END IF
  CALL StopTimer("FillFromTriplet")

  CALL DestructMatrix(local_matrix)
  CALL DestructTripletList(sorted_triplet_list)
