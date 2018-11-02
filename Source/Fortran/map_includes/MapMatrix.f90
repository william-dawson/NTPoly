  !! Convert to a triplet list, map the triplet list, fill.
  CALL ConstructEmptyMatrix(outmat, inmat)
  CALL GetMatrixTripletList(inmat, inlist)
  CALL MapTripletList(inlist, outlist, proc)
  CALL FillMatrixFromTripletList(outmat, inlist)

  !! Cleanup
  CALL DestructTripletList(inlist)
  CALL DestructTripletList(outlist)
