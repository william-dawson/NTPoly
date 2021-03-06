  !! Local Data
  INTEGER :: II, JJ
  INTEGER :: total

  !! There can't be more than one entry per row
  CALL ConstructTripletList(triplet_list, this%local_rows)

  total = 0
  !! Find local identity values
  DO JJ = this%start_row, this%end_row - 1
     DO II = this%start_column, this%end_column - 1
        IF (JJ .EQ. II .AND. JJ .LE. this%actual_matrix_dimension) THEN
           total = total + 1
           triplet_list%DATA(total)%index_column = II
           triplet_list%DATA(total)%index_row = JJ
           triplet_list%DATA(total)%point_value = 1.0
        END IF
     END DO
  END DO
  triplet_list%CurrentSize = total

  !! Finish constructing
  CALL FillMatrixFromTripletList(this, triplet_list, prepartitioned_in=.TRUE.)

  !! Cleanup
  CALL DestructTripletList(triplet_list)
