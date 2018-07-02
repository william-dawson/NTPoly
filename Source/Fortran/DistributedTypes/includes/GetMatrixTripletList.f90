  !! Merge all the local data
  CALL MergeMatrixLocalBlocks(this, merged_local_data)

  CALL MatrixToTripletList(merged_local_data, triplet_list)
  CALL ShiftTripletList(triplet_list, this%start_row - 1, &
       & this%start_column - 1)
