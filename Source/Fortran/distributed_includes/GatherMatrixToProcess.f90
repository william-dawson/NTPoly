  !! Local Variables
  INTEGER :: list_size
  INTEGER :: mat_dim
  INTEGER :: II

  !! Setup
  mat_dim = this%actual_matrix_dimension

  !! Local List
  CALL GetMatrixTripletList(this, tlist)
  list_size = tlist%CurrentSize

  !! Send this to the target
  ALLOCATE(slist(this%process_grid%slice_size))
  DO II = 1, this%process_grid%slice_size
     CALL ConstructTripletList(slist(II))
  END DO
  CALL ConstructTripletList(slist(within_slice_id + 1), list_size)
  slist(within_slice_id+1)%DATA(:list_size) = tlist%DATA(:list_size)
  CALL DestructTripletList(tlist)
  CALL RedistributeTripletLists(slist, this%process_grid%within_slice_comm, &
       & tlist)

  !! Create the local matrix
  IF (this%process_grid%within_slice_rank .EQ. within_slice_id) THEN
     CALL SortTripletList(tlist, mat_dim, mat_dim, sorted, .TRUE.)
     CALL ConstructMatrixFromTripletList(local_mat, sorted, &
          & mat_dim, mat_dim)
  END IF

  !! Cleanup
  CALL DestructTripletList(tlist)
  CALL DestructTripletList(sorted)
  DO II = 1, this%process_grid%slice_size
     CALL DestructTripletList(slist(II))
  END DO
