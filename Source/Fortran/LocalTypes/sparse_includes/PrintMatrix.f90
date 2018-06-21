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

  CALL this%PrintHeader(file_handler)

  WRITE(file_handler,'(A)') "%"
  WRITE(file_handler,*) this%rows, this%columns, size_of_this
  DO counter = 1,size_of_this
    CALL triplet_list%data(counter)%WriteToFile(file_handler)
  END DO

  IF (PRESENT(file_name_in)) THEN
     CLOSE(file_handler)
  END IF

  CALL triplet_list%Destruct()
