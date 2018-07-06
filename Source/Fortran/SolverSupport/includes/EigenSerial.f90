  !! Local Data
  INTEGER :: counter, list_size
  INTEGER :: mat_dim

  mat_dim = this%actual_matrix_dimension

  !! Gather on a single processor
  CALL GetMatrixTripletList(this, triplet_list)
  ALLOCATE(send_list(this%process_grid%slice_size))
  CALL ConstructTripletList(send_list(1), triplet_list%CurrentSize)
  DO counter = 2, this%process_grid%slice_size
     CALL ConstructTripletList(send_list(counter))
  END DO
  list_size = triplet_list%CurrentSize
  send_list(1)%data(:list_size) = triplet_list%data(:list_size)
  CALL DestructTripletList(triplet_list)
  CALL RedistributeTripletLists(send_list, &
       & this%process_grid%within_slice_comm, triplet_list)

  !! Perform the local decomposition
  CALL ConstructTripletList(triplet_w)
  IF (this%process_grid%within_slice_rank .EQ. 0) THEN
     CALL SortTripletList(triplet_list, mat_dim, mat_dim, &
          & sorted_triplet_list, .TRUE.)
     CALL ConstructMatrixFromTripletList(local_a, sorted_triplet_list, &
          & mat_dim, mat_dim)

     CALL ConstructMatrixDFromS(local_a, dense_a)
     IF (PRESENT(eigenvalues_out)) THEN
        CALL EigenDecomposition(dense_a, dense_v, dense_w)
        CALL ConstructTripletList(triplet_w, mat_dim)
        DO counter = 1, mat_dim
           triplet_w%data(counter)%index_row = counter
           triplet_w%data(counter)%index_column = counter
           triplet_w%data(counter)%point_value = dense_w%data(counter,1)
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
  CALL DestructTripletList(triplet_list)
  CALL DestructTripletList(triplet_w)
  CALL DestructTripletList(sorted_triplet_list)
  DO counter = 1, this%process_grid%slice_size
     CALL DestructTripletList(send_list(counter))
  END DO
  DEALLOCATE(send_list)
