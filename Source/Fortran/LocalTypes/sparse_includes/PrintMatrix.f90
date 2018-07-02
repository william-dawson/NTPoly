  INTEGER :: file_handler
  INTEGER :: counter
  INTEGER :: size_of_this

  !! Process Optional Parameters
  IF (PRESENT(file_name_in)) THEN
     file_handler = 16
     OPEN(unit = file_handler, file = file_name_in)
  ELSE
     file_handler = 6
  END IF

  !! Print
  CALL MatrixToTripletList(this,triplet_list)

  size_of_this = this%outer_index(this%columns+1)

#ifdef ISCOMPLEX
  WRITE(file_handler,'(A)') "%%MatrixMarket matrix coordinate complex general"
#else
  WRITE(file_handler,'(A)') "%%MatrixMarket matrix coordinate real general"
#endif

  WRITE(file_handler,'(A)') "%"
  WRITE(file_handler,*) this%rows, this%columns, size_of_this
  DO counter = 1,size_of_this
#ifdef ISCOMPLEX
     WRITE(file_handler,*) triplet_list%data(counter)%index_row, &
          & triplet_list%data(counter)%index_column, &
          & REAL(triplet_list%data(counter)%point_value), &
          & AIMAG(triplet_list%data(counter)%point_value)
#else
     WRITE(file_handler,*) triplet_list%data(counter)%index_row, &
          & triplet_list%data(counter)%index_column, &
          & triplet_list%data(counter)%point_value
#endif
  END DO

  IF (PRESENT(file_name_in)) THEN
     CLOSE(file_handler)
  END IF
  CALL DestructTripletList(triplet_list)
