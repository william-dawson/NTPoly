  !! Local data
  INTEGER :: working_col
  INTEGER :: values_per_hash
  INTEGER :: PII, HII, RII, CII
  INTEGER :: local_pii

  !! Loop over the hash structure — parallelize over RII (rows).
  !! Each RII reads from column RII of the 2D arrays. PII is a shared
  !! counter; use atomic capture to get a unique slot in pruned_list.
  PII = 1
!$acc parallel loop independent private(local_pii)
  DO RII = 1, mat_c_rows
     DO CII = 1, (mat_c_columns - 1) / memorypool%hash_size + 1
        !! Sort the elements in a hash
        values_per_hash = memorypool%inserted_per_bucket(CII, RII)
        memorypool%inserted_per_bucket(CII, RII) = 0
        !! Copy them
        DO HII = 1, values_per_hash
           working_col = memorypool%hash_index(HII + &
                & (CII - 1) * memorypool%hash_size, RII)
           working_value = &
                & memorypool%value_array(working_col, RII)
           memorypool%value_array(working_col, RII) = 0
           memorypool%dirty_array(working_col, RII) = .FALSE.
           !! If above threshold, insert
           IF (ABS(alpha*working_value) .GT. threshold) THEN
!$acc atomic capture
              local_pii = PII
              PII = PII + 1
!$acc end atomic
              memorypool%pruned_list(local_pii)%point_value = alpha*working_value
              memorypool%pruned_list(local_pii)%index_column = working_col
              memorypool%pruned_list(local_pii)%index_row = RII
           END IF
        END DO
     END DO
  END DO

  !! Convert to matrix (CPU — sort + CSR construction)
  CALL ConstructTripletList(unsorted_pruned_list, PII - 1)
  unsorted_pruned_list%DATA(:) = memorypool%pruned_list(1:PII - 1)
  CALL SortTripletList(unsorted_pruned_list, mat_c_columns, mat_c_rows, &
       & sorted_pruned_list, bubble_in = .TRUE.)
  CALL ConstructMatrixFromTripletList(matAB, sorted_pruned_list, mat_c_rows, &
       & mat_c_columns)
  CALL DestructTripletList(sorted_pruned_list)
  CALL DestructTripletList(unsorted_pruned_list)
