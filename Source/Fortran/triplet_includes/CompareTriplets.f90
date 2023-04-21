  IF (tripA%index_column .GT. tripB%index_column) THEN
     islessthan = .TRUE.
  ELSE IF (tripA%index_column .EQ. tripB%index_column .AND. &
       & tripA%index_row .GT. tripB%index_row) THEN
     islessthan = .TRUE.
  ELSE
     islessthan = .FALSE.
  END IF
