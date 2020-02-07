  CALL ComposeMatrix(this, merged)
  CALL MatrixToTripletList(merged, triplet_list)
  CALL DestructMatrix(merged)