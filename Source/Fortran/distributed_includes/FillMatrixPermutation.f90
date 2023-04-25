  !! Local Data
  INTEGER :: II, KK

  !! Build Local Triplet List
  !! There can't be more than one entry per row
  CALL ConstructTripletList(tlist, this%local_rows)
  KK = 0
  IF (rows) THEN
     DO II=this%start_row,this%end_row-1
        IF (permutation_vector(II) .GE. this%start_column .AND. &
             & permutation_vector(II) .LT. this%end_column) THEN
           KK = KK + 1
           tlist%DATA(KK)%index_column = permutation_vector(II)
           tlist%DATA(KK)%index_row = II
           tlist%DATA(KK)%point_value = 1.0
        END IF
     END DO
  ELSE
     DO II=this%start_column,this%end_column-1
        IF (permutation_vector(II) .GE. this%start_row .AND. &
             & permutation_vector(II) .LT. this%end_row) THEN
           KK = KK + 1
           tlist%DATA(KK)%index_column = II
           tlist%DATA(KK)%index_row = permutation_vector(II)
           tlist%DATA(KK)%point_value = 1.0
        END IF
     END DO
  END IF
  tlist%CurrentSize = KK

  !! Finish constructing
  CALL FillMatrixFromTripletList(this, tlist, prepartitioned_in = .TRUE.)

  !! Cleanup
  CALL DestructTripletList(tlist)
