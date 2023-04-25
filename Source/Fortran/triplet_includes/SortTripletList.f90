  !! Local Data
  LOGICAL :: bubble
  LOGICAL :: swap_occured
  INTEGER, DIMENSION(:), ALLOCATABLE :: values_per_row
  INTEGER, DIMENSION(:), ALLOCATABLE :: offset_array
  INTEGER, DIMENSION(:), ALLOCATABLE :: inserted_per_row
  !! Counters and temporary variables
  INTEGER :: II, idx
  INTEGER :: alloc_stat
  INTEGER :: list_length

  IF (PRESENT(bubble_in)) THEN
     bubble = bubble_in
  ELSE
     bubble = .TRUE.
  END IF

  list_length = input_list%CurrentSize

  IF (bubble .AND. list_length .GT. matrix_rows*matrix_columns * 0.1) THEN
     CALL SortDenseTripletList(input_list, matrix_columns, matrix_rows, &
          & sorted_list)
  ELSE
     !! Data Allocation
     CALL ConstructTripletList(sorted_list, list_length)
     ALLOCATE(values_per_row(matrix_columns), stat = alloc_stat)
     ALLOCATE(offset_array(matrix_columns), stat = alloc_stat)
     ALLOCATE(inserted_per_row(matrix_columns), stat = alloc_stat)

     !! Initial one dimensional sort
     values_per_row = 0
     inserted_per_row = 0

     !! Do a first pass bucket sort
     DO II = 1, input_list%CurrentSize
        values_per_row(input_list%DATA(II)%index_column) = &
             & values_per_row(input_list%DATA(II)%index_column) + 1
     END DO
     offset_array(1) = 1
     DO II = 2, UBOUND(offset_array, dim = 1)
        offset_array(II) = offset_array(II - 1) + &
             & values_per_row(II - 1)
     END DO
     DO II = 1, input_list%CurrentSize
        idx = input_list%DATA(II)%index_column
        sorted_list%DATA(offset_array(idx) + inserted_per_row(idx)) = &
             & input_list%DATA(II)
        inserted_per_row(idx) = inserted_per_row(idx) + 1
     END DO

     !! Finish with bubble sort
     !! Not necessary for transposing or unpacking.
     swap_occured = .TRUE.
     IF (bubble) THEN
        DO WHILE (swap_occured .EQV. .TRUE.)
           swap_occured = .FALSE.
           DO II = 2, sorted_list%CurrentSize
              IF (CompareTriplets(sorted_list%DATA(II - 1), &
                   & sorted_list%DATA(II))) THEN
                 trip = sorted_list%DATA(II)
                 sorted_list%DATA(II) = sorted_list%DATA(II - 1)
                 sorted_list%DATA(II - 1) = trip
                 swap_occured = .TRUE.
              END IF
           END DO
        END DO
     END IF

     !! Cleanup
     DEALLOCATE(values_per_row)
     DEALLOCATE(offset_array)
     DEALLOCATE(inserted_per_row)
  END IF
