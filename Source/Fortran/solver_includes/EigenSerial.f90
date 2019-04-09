  !! Local Variables
  INTEGER :: II, mat_dim

  mat_dim = this%actual_matrix_dimension
  CALL CopyMatrix(this, tempmat)
  CALL GatherMatrixToProcess(tempmat, local_a, 0)

  !! Perform the local decomposition
  CALL ConstructTripletList(triplet_w)
  CALL ConstructTripletList(triplet_list)
  IF (this%process_grid%within_slice_rank .EQ. 0) THEN
     CALL ConstructMatrixDFromS(local_a, dense_a)
     IF (PRESENT(eigenvalues_out)) THEN
        CALL EigenDecomposition(dense_a, dense_v, dense_w)
        CALL ConstructTripletList(triplet_w, mat_dim)
        DO II = 1, mat_dim
           triplet_w%data(II)%index_row = II
           triplet_w%data(II)%index_column = II
           triplet_w%data(II)%point_value = dense_w%data(II,1)
        END DO
     ELSE
        CALL EigenDecomposition(dense_a, dense_v)
     END IF

     CALL ConstructMatrixSFromD(dense_v, local_v, fixed_params%threshold)
     CALL MatrixToTripletList(local_v, triplet_list)
  END IF

  !! Build The Full Matrices
  CALL ConstructEmptyMatrix(eigenvectors, this)
  CALL FillMatrixFromTripletList(eigenvectors, triplet_list, .TRUE.)

  IF (PRESENT(eigenvalues_out)) THEN
     CALL ConstructEmptyMatrix(eigenvalues_out, this)
     CALL FillMatrixFromTripletList(eigenvalues_out, triplet_w, .TRUE.)
  END IF

  !! Cleanup
  CALL DestructMatrix(dense_a)
  CALL DestructMatrix(dense_v)
  CALL DestructMatrix(dense_w)
  CALL DestructMatrix(sparse)
  CALL DestructTripletList(triplet_w)
  CALL DestructTripletList(triplet_list)
  CALL DestructMatrix(tempmat)
