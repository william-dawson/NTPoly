
  INTEGER :: II
  LOGICAL :: valid
  INTEGER :: num_slices
  INTEGER :: my_slice

  IF (PRESENT(num_slices_in)) THEN
     num_slices = num_slices_in
     IF (PRESENT(my_slice_in)) THEN
        my_slice = my_slice_in
     ELSE
        my_slice = 0
     END IF
  ELSE
     num_slices = 1
     my_slice = 1
  END IF

  CALL ConstructTripletList(outlist)
  DO II = my_slice+1, inlist%CurrentSize, num_slices
     CALL GetTripletAt(inlist, II, temp)
     CALL proc(temp%index_row, temp%index_column, temp%point_value, valid)
     IF (valid) THEN
        CALL AppendToTripletList(outlist, temp)
     END IF
  END DO
