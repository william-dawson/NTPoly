  INTEGER :: inserted_vals
  INTEGER :: idx_a, idx_b, idx_hash
  INTEGER :: elements_per_inner_a, elements_per_inner_b
  LOGICAL :: is_set
  !! Counters
  INTEGER :: II, AA, BB

  !! Multiply
  DO II = 1, matAT%columns
     elements_per_inner_a = matAT%outer_index(II + 1) - &
          & matAT%outer_index(II)
     DO AA = 1, elements_per_inner_a
        val_a = matAT%values(matAT%outer_index(II) + AA)
        idx_a = matAT%inner_index(matAT%outer_index(II) + AA)
        elements_per_inner_b = matBT%outer_index(idx_a + 1) - &
             & matBT%outer_index(idx_a)
        DO BB = 1, elements_per_inner_b
           idx_b = matBT%inner_index(matBT%outer_index(idx_a) + BB)
           val_b = matBT%values(matBT%outer_index(idx_a)+ BB)
           val_c = memorypool%value_array(idx_b, II)
           is_set = memorypool%dirty_array(idx_b, II)
           IF (is_set .EQV. .FALSE.) THEN
              memorypool%dirty_array(idx_b, II) = .TRUE.
              idx_hash = (idx_b - 1) / memorypool%hash_size
              inserted_vals = & 
                   & memorypool%inserted_per_bucket(idx_hash + 1, II) + 1
              memorypool%inserted_per_bucket(idx_hash + 1, II) = &
                   & inserted_vals
              memorypool%hash_index(&
                   & inserted_vals + idx_hash * memorypool%hash_size, &
                   & II) = idx_b
           END IF
           memorypool%value_array(idx_b, II) = val_c + val_a * val_b
        END DO
     END DO
  END DO
