  !! Local Variables
  TYPE(TLISTTYPE) :: tlist, slist
  TYPE(TTYPE) :: triplet
  !! Temporary Variables
  INTEGER :: II
  INTEGER :: new_dim

  !! Get a triplet list with the values
  CALL GetMatrixTripletList(this, tlist)
  CALL ConstructTripletList(slist)

  !! Filter and shift the triplet list
  DO II = 1, tlist%CurrentSize
     CALL GetTripletAt(tlist, II, triplet)
     IF (triplet%index_row .GE. start_row .AND. &
          & triplet%index_row .LE. end_row .AND. &
          & triplet%index_column .GE. start_column .AND. &
          & triplet%index_column .LE. end_column) THEN
        triplet%index_row = triplet%index_row - start_row + 1
        triplet%index_column = triplet%index_column - start_column + 1
        CALL AppendToTripletList(slist, triplet)
     END IF
  END DO

  new_dim = MAX(end_row - start_row + 1, end_column - start_column + 1)
  CALL ConstructEmptyMatrix(submatrix, new_dim, &
       & process_grid_in = this%process_grid, is_complex_in = this%is_complex)
  CALL FillMatrixFromTripletList(submatrix, slist, preduplicated_in = .TRUE.)

  !! Cleanup
  CALL DestructTripletList(tlist)
  CALL DestructTripletList(slist)
