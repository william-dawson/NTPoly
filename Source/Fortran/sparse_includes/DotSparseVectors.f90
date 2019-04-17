  INTEGER :: working_index_a, working_index_b
  !! Counter Variables
  INTEGER :: counter_a, counter_b
  INTEGER :: sizea, sizeb

  counter_a = 1
  counter_b = 1
  sizea = SIZE(inner_index_a)
  sizeb = SIZE(inner_index_b)
  product = 0
  sum_loop: DO WHILE(counter_a .LE. sizea .AND. counter_b .LE. sizeb)
     !! Current inner indices and values
     working_index_a = inner_index_a(counter_a)
     working_index_b = inner_index_b(counter_b)
     !! Figure out from which vector an insertion will be performed
     IF (working_index_a .EQ. working_index_b) THEN
        working_value_a = values_a(counter_a)
        working_value_b = values_b(counter_b)
        product = product + working_value_a * working_value_b
        counter_a = counter_a + 1
        counter_b = counter_b + 1
     ELSE IF (working_index_a .GT. working_index_b) THEN
        DO WHILE(counter_b .LE. sizeb)
          counter_b = counter_b + 1
          IF (inner_index_b(counter_b) .GE. working_index_a) EXIT
        END DO
     ELSE !! implies working_index_b > working_index_b
        DO WHILE(counter_a .LE. sizea)
          counter_a = counter_a + 1
          IF (inner_index_a(counter_a) .GE. working_index_b) EXIT
        END DO
     END IF
  END DO sum_loop
