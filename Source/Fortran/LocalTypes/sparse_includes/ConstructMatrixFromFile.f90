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
  CALL triplet_list%Init

  !! Read Values
  DO counter = 1, temp_total_values
     CALL temporary%ReadFromFile(file_handler)
     CALL triplet_list%Append(temporary)
  END DO

  CLOSE(file_handler)
  CALL triplet_list%Symmetrize(pattern_type)
  CALL SortTripletList(triplet_list, temp_columns, temp_rows, &
       & sorted_triplet_list)
  CALL ConstructMatrixFromTripletList(this, sorted_triplet_list, temp_rows, &
       & temp_columns)

  CALL triplet_list%Destruct
  CALL sorted_triplet_list%Destruct
