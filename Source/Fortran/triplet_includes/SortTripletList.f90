  !! Local Data
  LOGICAL :: bubble
  INTEGER, DIMENSION(:), ALLOCATABLE :: values_per_row
  INTEGER, DIMENSION(:), ALLOCATABLE :: offset_array
  INTEGER, DIMENSION(:), ALLOCATABLE :: inserted_per_row
  !! Counters and temporary variables
  INTEGER :: II, JJ, idx
  INTEGER :: alloc_stat
  INTEGER :: list_length
  !! Per-bucket insertion sort variables
  INTEGER :: bucket_start, bucket_end, bucket_size

  IF (PRESENT(bubble_in)) THEN
     bubble = bubble_in
  ELSE
     bubble = .TRUE.
  END IF

  list_length = input_list%CurrentSize

  IF (bubble .AND. REAL(list_length) .GT. &
       & REAL(matrix_rows)*REAL(matrix_columns) * 0.1) THEN
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

     !! Sort within each column bucket by row.
     !! The old code used a bubble sort over the ENTIRE list
     !! (O(max_bucket * nz)) which thrashes L3 cache because it scans
     !! millions of elements up to ~76 times.  Instead, do an insertion
     !! sort per bucket — each bucket is small (~75 entries) and fits
     !! entirely in L1 cache, making the sort ~50x faster.
     IF (bubble) THEN
        DO idx = 1, matrix_columns
           bucket_start = offset_array(idx)
           bucket_size = values_per_row(idx)
           bucket_end = bucket_start + bucket_size - 1
           IF (bucket_size .GT. 1) THEN
              DO II = bucket_start + 1, bucket_end
                 trip = sorted_list%DATA(II)
                 JJ = II - 1
                 DO WHILE (JJ .GE. bucket_start)
                    IF (sorted_list%DATA(JJ)%index_row <= trip%index_row) EXIT
                    sorted_list%DATA(JJ + 1) = sorted_list%DATA(JJ)
                    JJ = JJ - 1
                 END DO
                 sorted_list%DATA(JJ + 1) = trip
              END DO
           END IF
        END DO
     END IF

     !! Cleanup
     DEALLOCATE(values_per_row)
     DEALLOCATE(offset_array)
     DEALLOCATE(inserted_per_row)
  END IF
