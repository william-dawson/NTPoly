  !! Gather as dense
  CALL GatherMatrixToProcess(this, local_s)
  CALL ConstructMatrixDFromS(local_s, local_d)

  !! Decompose
  CALL EigenDecomposition(local_d, V, W)

  !! Filter if necessary
  IF (nvals+1 .LE. V%rows) THEN
     V%DATA(:, nvals+1:) = 0
     W%DATA(nvals+1:, :) = 0
     W%DATA(:, nvals+1:) = 0
  END IF

  !! Convert results to triplet lists
  CALL ConstructMatrixSFromD(V, V_s, threshold_in=threshold)
  CALL ConstructMatrixSFromD(W, W_s, threshold_in=threshold)
  CALL MatrixToTripletList(V_s, V_t)
  CALL MatrixToTripletList(W_s, W_t)

  !! Distribute
  CALL ConstructEmptyMatrix(eigenvalues, this)
  IF (eigenvalues%process_grid%within_slice_rank .NE. 0) THEN
     CALL ConstructTripletList(W_t)
  END IF
  CALL FillMatrixFromTripletList(eigenvalues, W_t, preduplicated_in=.TRUE.)

  IF (PRESENT(eigenvectors_in)) THEN
     CALL ConstructEmptyMatrix(eigenvectors_in, this)
     IF (eigenvectors_in%process_grid%within_slice_rank .NE. 0) THEN
        CALL ConstructTripletList(V_t)
     END IF
     CALL FillMatrixFromTripletList(eigenvectors_in, V_t, &
          & preduplicated_in=.TRUE.)
  END IF

  !! Cleanup
  CALL DestructMatrix(local_s)
  CALL DestructMatrix(local_d)
  CALL DestructTripletList(V_t)
  CALL DestructTripletList(W_t)
