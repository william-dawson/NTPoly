  !! Local Data
  REAL(NTREAL) :: alpha
  REAL(NTREAL) :: threshold
  !! Temporary Variables
  INTEGER :: working_index_a, working_index_b
  !! Counter Variables
  INTEGER :: AA, BB, CC

  !! Process Optional Parameters
  IF (.NOT. PRESENT(alpha_in)) THEN
     alpha = 1.0_NTREAL
  ELSE
     alpha = alpha_in
  END IF
  IF (.NOT. PRESENT(threshold_in)) THEN
     threshold = 0.0_NTREAL
  ELSE
     threshold = threshold_in
  END IF

  AA = 1
  BB = 1
  CC = 1
  DO WHILE(AA .LE. SIZE(inner_index_a) .AND. BB .LE. SIZE(inner_index_b))
     !! Current inner indices and values
     working_index_a = inner_index_a(AA)
     working_index_b = inner_index_b(BB)
     working_value_a = alpha*values_a(AA)
     working_value_b = values_b(BB)
     !! Figure out from which vector an insertion will be performed
     IF (working_index_a .EQ. working_index_b) THEN
        IF (ABS(working_value_a + working_value_b) .GT. threshold) THEN
           inner_index_c(CC) = working_index_a
           values_c(CC) = working_value_a + working_value_b
           CC = CC + 1
        END IF
        AA = AA + 1
        BB = BB + 1
     ELSE IF (working_index_a .GT. working_index_b) THEN
        IF (ABS(working_value_b) .GT. threshold) THEN
           inner_index_c(CC) = working_index_b
           values_c(CC) = working_value_b
           CC = CC + 1
        END IF
        BB = BB + 1
     ELSE !! implies working_index_b > working_index_a
        IF (ABS(working_value_a) .GT. threshold) THEN
           inner_index_c(CC) = working_index_a
           values_c(CC) = working_value_a
           CC = CC + 1
        END IF
        AA = AA + 1
     END IF
  END DO

  !! Handle case where one was blank
  DO WHILE (AA .LE. SIZE(inner_index_a))
     inner_index_c(CC) = inner_index_a(AA)
     values_c(CC) = values_a(AA)*alpha
     AA = AA + 1
     CC = CC + 1
  END DO
  DO WHILE (BB .LE. SIZE(inner_index_b))
     inner_index_c(CC) = inner_index_b(BB)
     values_c(CC) = values_b(BB)
     BB = BB + 1
     CC = CC + 1
  END DO

  total_values_c = CC - 1
