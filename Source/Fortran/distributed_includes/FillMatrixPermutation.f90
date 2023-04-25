  !! Local Data
  INTEGER :: II, KK

  !! Build Local Triplet List
  !! There can't be more than one entry per row
  CALL ConstructTripletList(triplet_list, this%local_rows)
  KK = 0
  IF (rows) THEN
     DO II=this%start_row,this%end_row-1
        IF (permutation_vector(II) .GE. this%start_column .AND. &
             & permutation_vector(II) .LT. this%end_column) THEN
           KK = KK + 1
           triplet_list%DATA(KK)%index_column = permutation_vector(II)
           triplet_list%DATA(KK)%index_row = II
           triplet_list%DATA(KK)%point_value = 1.0
        END IF
     END DO
  ELSE
     DO II=this%start_column,this%end_column-1
        IF (permutation_vector(II) .GE. this%start_row .AND. &
             & permutation_vector(II) .LT. this%end_row) THEN
           KK = KK + 1
           triplet_list%DATA(KK)%index_column = II
           triplet_list%DATA(KK)%index_row = permutation_vector(II)
           triplet_list%DATA(KK)%point_value = 1.0
        END IF
     END DO
  END IF
  triplet_list%CurrentSize = KK

  !! Finish constructing
  CALL FillMatrixFromTripletList(this, triplet_list, prepartitioned_in=.TRUE.)

  !! Cleanup
  CALL DestructTripletList(triplet_list)
