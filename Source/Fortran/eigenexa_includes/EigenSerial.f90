!! Gather as dense
CALL GatherMatrixToProcess(this, local_s)
CALL ConstructMatrixDFromS(local_s, local_d)

!! Decompose
CALL EigenDecomposition(local_d, V, W)

WRITE(*,*) ":::::", W%DATA

!! Convert results to triplet lists
CALL ConstructMatrixSFromD(V, V_s, threshold_in=threshold)
CALL ConstructMatrixSFromD(W, W_s, threshold_in=threshold)
CALL MatrixToTripletList(V_s, V_t)
CALL MatrixToTripletList(W_s, W_t)

!! Distribute
CALL ConstructEmptyMatrix(eigenvectors, this)
CALL ConstructEmptyMatrix(eigenvalues, this)
IF (eigenvectors%process_grid%within_slice_rank .NE. 0) THEN
   CALL ConstructTripletList(V_t)
   CALL ConstructTripletList(W_t)
END IF
CALL FillMatrixFromTripletList(eigenvectors, V_t, preduplicated_in=.TRUE.)
CALL FillMatrixFromTripletList(eigenvalues, W_t, preduplicated_in=.TRUE.)

CALL PrintMatrix(eigenvalues)

!! Cleanup
CALL DestructMatrix(local_s)
CALL DestructMatrix(local_d)
CALL DestructTripletList(V_t)
CALL DestructTripletList(W_t)