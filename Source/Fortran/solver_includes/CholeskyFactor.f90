  !! Convert to dense and factorize
  CALL CopyMatrix(mat, tempmat)
  CALL GatherMatrixToProcess(tempmat, smat, 0)
  CALL ConstructTripletList(tlist)
  IF (mat%process_grid%within_slice_rank .EQ. 0) THEN
     CALL ConstructMatrixDFromS(smat, dmat)
     CALL CholeskyFactor(dmat, dfac)
     CALL ConstructMatrixSFromD(dfac, sfac, parameters%threshold)
     CALL MatrixToTripletList(sfac, tlist)
  END IF

  !! Build The Full Matrices
  CALL ConstructEmptyMatrix(factor, mat)
  CALL FillMatrixFromTripletList(factor, tlist, .TRUE.)

  !! Cleanup
  CALL DestructMatrix(smat)
  CALL DestructMatrix(dmat)
  CALL DestructMatrix(sfac)
  CALL DestructMatrix(dfac)
  CALL DestructMatrix(tempmat)
