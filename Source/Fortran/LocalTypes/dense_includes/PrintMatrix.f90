  !! Local Data
  INTEGER :: file_handler
  INTEGER II, JJ

  !! Process Optional Parameters
  IF (PRESENT(file_name_in)) THEN
     file_handler = 16
     OPEN(unit = file_handler, file = file_name_in)
  ELSE
     file_handler = 6
  END IF

  CALL this%PrintHeader(file_handler)
  DO II = 1, this%columns
     DO JJ = 1, this%rows
        WRITE(file_handler,*) this%data(II,JJ)
     END DO
  END DO

  IF (PRESENT(file_name_in)) THEN
     CLOSE(file_handler)
  END IF
