  !! Local Data
  LOGICAL :: bubble
  LOGICAL :: swap_occured
  INTEGER, DIMENSION(:), ALLOCATABLE :: values_per_row
  INTEGER, DIMENSION(:), ALLOCATABLE :: offset_array
  INTEGER, DIMENSION(:), ALLOCATABLE :: inserted_per_row
  !! Counters and temporary variables
  INTEGER :: counter
  INTEGER :: temp_index
  INTEGER :: alloc_stat
  INTEGER :: list_length

  IF (PRESENT(bubble_in)) THEN
     bubble = bubble_in
  ELSE
     bubble = .TRUE.
  END IF

  list_length = input_list%CurrentSize

  IF (bubble .AND. list_length .GT. matrix_rows*matrix_columns*0.1) THEN
     CALL SortDenseTripletList(input_list, matrix_columns, matrix_rows, &
          & sorted_list)
  ELSE
     !! Data Allocation
     sorted_list = ConstructTripletList(list_length)
     ALLOCATE(values_per_row(matrix_columns), stat=alloc_stat)
     ALLOCATE(offset_array(matrix_columns), stat=alloc_stat)
     ALLOCATE(inserted_per_row(matrix_columns), stat=alloc_stat)

     !! Initial one dimensional sort
     values_per_row = 0
     inserted_per_row = 0

     !! Do a first pass bucket sort
     DO counter = 1, input_list%CurrentSize
        values_per_row(input_list%data(counter)%index_column) = &
             & values_per_row(input_list%data(counter)%index_column) + 1
     END DO
     offset_array(1) = 1
     DO counter = 2, UBOUND(offset_array,dim=1)
        offset_array(counter) = offset_array(counter-1) + &
             & values_per_row(counter-1)
     END DO
     DO counter = 1, input_list%CurrentSize
        temp_index = input_list%data(counter)%index_column
        sorted_list%data(offset_array(temp_index)+inserted_per_row(temp_index))=&
             & input_list%data(counter)
        inserted_per_row(temp_index) = inserted_per_row(temp_index) + 1
     END DO

     !! Finish with bubble sort
     !! Not necessary for transposing or unpacking.
     swap_occured = .TRUE.
     IF (bubble) THEN
        DO WHILE (swap_occured .EQV. .TRUE.)
           swap_occured = .FALSE.
           DO counter = 2, sorted_list%CurrentSize
              IF (CompareTriplets(sorted_list%data(counter-1), &
                   & sorted_list%data(counter))) THEN
                 temporary = sorted_list%data(counter)
                 sorted_list%data(counter) = sorted_list%data(counter-1)
                 sorted_list%data(counter-1) = temporary
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
