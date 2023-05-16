  INTEGER :: working_index_a, working_index_b
  !! Counter Variables
  INTEGER :: AA, BB, CC

  AA = 1
  BB = 1
  CC = 1
  DO WHILE(AA .LE. SIZE(inner_index_a) .AND. BB .LE. SIZE(inner_index_b))
     !! Current inner indices and values
     working_index_a = inner_index_a(AA)
     working_index_b = inner_index_b(BB)
     working_value_a = values_a(AA)
     working_value_b = values_b(BB)
     !! Figure out from which vector an insertion will be performed
     IF (working_index_a .EQ. working_index_b) THEN
        inner_index_c(CC) = working_index_a
        values_c(CC) = working_value_a * working_value_b
        AA = AA + 1
        BB = BB + 1
        CC = CC + 1
     ELSE IF (working_index_a .GT. working_index_b) THEN
        BB = BB + 1
     ELSE !! implies working_index_b > working_index_b
        AA = AA + 1
     END IF
  END DO
  total_values_c = CC - 1
