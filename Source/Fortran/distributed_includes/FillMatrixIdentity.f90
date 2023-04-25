  !! Local Data
  INTEGER :: II, JJ, KK

  !! There can't be more than one entry per row
  CALL ConstructTripletList(triplet_list, this%local_rows)

  KK = 0
  !! Find local identity values
  DO JJ = this%start_row, this%end_row - 1
     DO II = this%start_column, this%end_column - 1
        IF (JJ .EQ. II .AND. JJ .LE. this%actual_matrix_dimension) THEN
           KK = KK + 1
           triplet_list%DATA(KK)%index_column = II
           triplet_list%DATA(KK)%index_row = JJ
           triplet_list%DATA(KK)%point_value = 1.0
        END IF
     END DO
  END DO
  triplet_list%CurrentSize = KK

  !! Finish constructing
  CALL FillMatrixFromTripletList(this, triplet_list, prepartitioned_in=.TRUE.)

  !! Cleanup
  CALL DestructTripletList(triplet_list)
