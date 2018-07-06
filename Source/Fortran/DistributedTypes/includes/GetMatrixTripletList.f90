  !! Merge all the local data
  CALL MergeMatrixLocalBlocks(working_matrix, merged_local_data)

  CALL MatrixToTripletList(merged_local_data, triplet_list)
  CALL ShiftTripletList(triplet_list, working_matrix%start_row - 1, &
       & working_matrix%start_column - 1)

  CALL DestructMatrix(working_matrix)
