  !! Local Data
  INTEGER :: II, KK, first_diagonal, last_diagonal

  !! There can not be more than one entry per row
  first_diagonal = MAX(this%start_row, this%start_column)
  last_diagonal = MIN(this%end_row - 1, this%end_column - 1, &
       & this%actual_matrix_dimension)
  CALL ConstructTripletList(tlist, MAX(0, last_diagonal - first_diagonal + 1))

  KK = 0
  !! Find local identity values
  DO II = first_diagonal, last_diagonal
     KK = KK + 1
     tlist%DATA(KK)%index_column = II
     tlist%DATA(KK)%index_row = II
     tlist%DATA(KK)%point_value = 1.0
  END DO
  tlist%CurrentSize = KK

  !! Finish constructing
  CALL FillMatrixFromTripletList(this, tlist, prepartitioned_in = .TRUE.)

  !! Cleanup
  CALL DestructTripletList(tlist)
