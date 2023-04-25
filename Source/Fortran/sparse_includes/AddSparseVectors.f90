  !! Local Data
  REAL(NTREAL) :: alpha
  REAL(NTREAL) :: threshold
  !! Temporary Variables
  INTEGER :: working_index_a, working_index_b
  !! Counter Variables
  INTEGER :: counter_a, counter_b, counter_c

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

  counter_a = 1
  counter_b = 1
  counter_c = 1
  sum_loop: DO WHILE(counter_a .LE. SIZE(inner_index_a) .AND. counter_b .LE. &
       & SIZE(inner_index_b))
     !! Current inner indices and values
     working_index_a = inner_index_a(counter_a)
     working_index_b = inner_index_b(counter_b)
     working_value_a = alpha*values_a(counter_a)
     working_value_b = values_b(counter_b)
     !! Figure out from which vector an insertion will be performed
     IF (working_index_a .EQ. working_index_b) THEN
        IF (ABS(working_value_a + working_value_b) .GT. threshold) THEN
           inner_index_c(counter_c) = working_index_a
           values_c(counter_c) = working_value_a + working_value_b
           counter_c = counter_c + 1
        END IF
        counter_a = counter_a + 1
        counter_b = counter_b + 1
     ELSE IF (working_index_a .GT. working_index_b) THEN
        IF (ABS(working_value_b) .GT. threshold) THEN
           inner_index_c(counter_c) = working_index_b
           values_c(counter_c) = working_value_b
           counter_c = counter_c + 1
        END IF
        counter_b = counter_b + 1
     ELSE !! implies working_index_b > working_index_a
        IF (ABS(working_value_a) .GT. threshold) THEN
           inner_index_c(counter_c) = working_index_a
           values_c(counter_c) = working_value_a
           counter_c = counter_c + 1
        END IF
        counter_a = counter_a + 1
     END IF
  END DO sum_loop

  !! Handle case where one was blank
  cleanup_a: DO WHILE (counter_a .LE. SIZE(inner_index_a))
     inner_index_c(counter_c) = inner_index_a(counter_a)
     values_c(counter_c) = values_a(counter_a)*alpha
     counter_a = counter_a + 1
     counter_c = counter_c + 1
  END DO cleanup_a
  cleanup_b: DO WHILE (counter_b .LE. SIZE(inner_index_b))
     inner_index_c(counter_c) = inner_index_b(counter_b)
     values_c(counter_c) = values_b(counter_b)
     counter_b = counter_b + 1
     counter_c = counter_c + 1
  END DO cleanup_b

  total_values_c = counter_c - 1
