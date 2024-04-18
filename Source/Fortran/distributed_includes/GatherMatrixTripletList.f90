  CALL GatherMatrixToProcess(this, lmat)
  CALL MatrixToTripletList(lmat, tlist)
  CALL DestructMatrix(lmat)
