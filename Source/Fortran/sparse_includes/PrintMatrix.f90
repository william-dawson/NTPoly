  !! Local Data
  INTEGER :: file_handler
  INTEGER :: counter
  INTEGER :: size_of_this
  CHARACTER(LEN=MAX_LINE_LENGTH) :: tempstr

  !! Process Optional Parameters
  IF (PRESENT(file_name_in)) THEN
     file_handler = 16
     OPEN(unit = file_handler, file = file_name_in)
  ELSE
     file_handler = 6
  END IF

  !! Print
  CALL MatrixToTripletList(this, triplet_list)

  size_of_this = this%outer_index(this%columns+1)

#ifdef ISCOMPLEX
  WRITE(file_handler,'(A)') "%%MatrixMarket matrix coordinate complex general"
#else
  WRITE(file_handler,'(A)') "%%MatrixMarket matrix coordinate real general"
#endif

  WRITE(file_handler,'(A)') "%"
  CALL WriteMMSize(tempstr, this%rows, this%columns, &
       & INT(size_of_this, KIND=NTLONG))
  WRITE(file_handler,'(A)') ADJUSTL(TRIM(tempstr))
  DO counter = 1,size_of_this
#ifdef ISCOMPLEX
     CALL WriteMMLine(tempstr, triplet_list%DATA(counter)%index_row, &
          & triplet_list%DATA(counter)%index_column, &
          & REAL(triplet_list%DATA(counter)%point_value), &
          & AIMAG(triplet_list%DATA(counter)%point_value))
     WRITE(file_handler,'(A)') ADJUSTL(TRIM(tempstr))
#else
     CALL WriteMMLine(tempstr, triplet_list%DATA(counter)%index_row, &
          & triplet_list%DATA(counter)%index_column, &
          & triplet_list%DATA(counter)%point_value)
     WRITE(file_handler,'(A)') ADJUSTL(TRIM(tempstr))
#endif
  END DO

  IF (PRESENT(file_name_in)) THEN
     CLOSE(file_handler)
  END IF
  CALL DestructTripletList(triplet_list)
