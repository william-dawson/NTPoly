  INTEGER :: working_index_a, working_index_b
  !! Counter Variables
  INTEGER :: counter_a, counter_b

  counter_a = 1
  counter_b = 1
  product = 0
  sum_loop: DO WHILE(counter_a .LE. SIZE(inner_index_a) .AND. counter_b .LE. &
       & SIZE(inner_index_b))
     !! Current inner indices and values
     working_index_a = inner_index_a(counter_a)
     working_index_b = inner_index_b(counter_b)
     working_value_a = values_a(counter_a)
     working_value_b = values_b(counter_b)
     !! Figure out from which vector an insertion will be performed
     IF (working_index_a .EQ. working_index_b) THEN
        product = product + working_value_a * working_value_b
        counter_a = counter_a + 1
        counter_b = counter_b + 1
     ELSE IF (working_index_a .GT. working_index_b) THEN
        counter_b = counter_b + 1
     ELSE !! implies working_index_b > working_index_b
        counter_a = counter_a + 1
     END IF
  END DO sum_loop
