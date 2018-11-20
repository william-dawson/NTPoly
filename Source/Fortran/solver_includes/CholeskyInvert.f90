  !! Convert to dense and factorize
  CALL GatherMatrixToProcess(mat, smat, 0)
  CALL ConstructTripletList(tlist)
  IF (mat%process_grid%within_slice_rank .EQ. 0) THEN
     CALL ConstructMatrixDFromS(smat, dmat)
     CALL CholeskYInvert(dmat, dinv)
     CALL ConstructMatrixSFromD(dinv, sinv, parameters%threshold)
     CALL MatrixToTripletList(sinv, tlist)
  END IF

  !! Build The Full Matrices
  CALL ConstructEmptyMatrix(factor, mat)
  CALL FillMatrixFromTripletList(factor, tlist, .TRUE.)

  !! Cleanup
  CALL DestructMatrix(smat)
  CALL DestructMatrix(dmat)
  CALL DestructMatrix(sinv)
  CALL DestructMatrix(dinv)
