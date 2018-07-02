  !! Local Data
  INTEGER :: temp_rows, temp_columns, temp_total_values
  CHARACTER(len=81) :: input_buffer
  INTEGER :: file_handler
  INTEGER :: counter
  LOGICAL :: found_comment_line
  LOGICAL :: error_occured
  file_handler = 16

  OPEN(file_handler,file=file_name,status='old')

  !! Parse the header.
  READ(file_handler,fmt='(A)') input_buffer
  error_occured = ParseMMHeader(input_buffer, sparsity_type, data_type, &
       & pattern_type)

  !! Extra Comment Lines
  found_comment_line = .TRUE.
  DO WHILE(found_comment_line)
     !read(file_handler,*) input_buffer
     READ(file_handler,fmt='(A)') input_buffer
     IF (.NOT. input_buffer(1:1) .EQ. '%') THEN
        found_comment_line = .FALSE.
     END IF
  END DO

  !! Main data
  READ(input_buffer,*) temp_rows, temp_columns, temp_total_values
  CALL ConstructTripletList(triplet_list)

  !! Read Values
  DO counter = 1, temp_total_values
#ifdef ISCOMPLEX
     READ(file_handler,*) temporary%index_row, temporary%index_column, &
          & real_val, comp_val
     temporary%point_value = CMPLX(real_val, comp_val, KIND=NTCOMPLEX)
#else
     READ(file_handler,*) temporary%index_row, temporary%index_column, &
          & temporary%point_value
#endif
     CALL AppendToTripletList(triplet_list,temporary)
  END DO

  CLOSE(file_handler)
  CALL SymmetrizeTripletList(triplet_list, pattern_type)
  CALL SortTripletList(triplet_list, temp_columns, temp_rows, &
       & sorted_triplet_list)
  CALL ConstructMatrixFromTripletList(this, sorted_triplet_list, temp_rows, &
       & temp_columns)

  CALL DestructTripletList(triplet_list)
  CALL DestructTripletList(sorted_triplet_list)
