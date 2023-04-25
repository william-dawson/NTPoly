  !! Local Data
  INTEGER :: rows, columns, total_values
  CHARACTER(len=81) :: input_buffer
  INTEGER, PARAMETER :: file_handler = 16
  LOGICAL :: found_comment_line
  LOGICAL :: error_occured
  INTEGER :: II

  OPEN(file_handler, file=file_name, status='old')

  !! Parse the header.
  READ(file_handler,fmt='(A)') input_buffer
  error_occured = ParseMMHeader(input_buffer, sparsity_type, data_type, &
       & pattern_type)

  !! Extra Comment Lines
  found_comment_line = .TRUE.
  DO WHILE(found_comment_line)
     READ(file_handler, fmt = '(A)') input_buffer
     IF (.NOT. input_buffer(1:1) .EQ. '%') THEN
        found_comment_line = .FALSE.
     END IF
  END DO

  !! Main data
  READ(input_buffer,*) rows, columns, total_values
  CALL ConstructTripletList(tlist)

  !! Read Values
  DO II = 1, total_values
#ifdef ISCOMPLEX
     READ(file_handler,*) temp%index_row, temp%index_column, real_val, comp_val
     temp%point_value = CMPLX(real_val, comp_val, KIND=NTCOMPLEX)
#else
     READ(file_handler,*) temp%index_row, temp%index_column, temp%point_value
#endif
     CALL AppendToTripletList(tlist, temp)
  END DO

  CLOSE(file_handler)
  CALL SymmetrizeTripletList(tlist, pattern_type)
  CALL SortTripletList(tlist, columns, rows, sorted_tlist)
  CALL ConstructMatrixFromTripletList(this, sorted_tlist, rows, columns)

  CALL DestructTripletList(tlist)
  CALL DestructTripletList(sorted_tlist)
