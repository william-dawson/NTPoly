  !! Local Variables
  TYPE(Matrix_ps) :: invecT
  TYPE(Matrix_ps) :: S
  TYPE(Matrix_ps) :: rotation
  TYPE(Matrix_ldr) :: local_w
  REAL(NTREAL) :: temp_val
  INTEGER :: II

  !! Compute the overlap matrix.
  CALL TransposeMatrix(invec, invecT)
  IF (invec%is_complex) THEN
     CALL ConjugateMatrix(invecT)
  END IF
  CALL MatrixMultiply(invecT, invec, S, memory_pool_in=pool)

  !! Slice out the overlap matrix.
  CALL ResizeMatrix(S, num_vecs)
  CALL GatherMatrixToProcess(S, local_sparse, 0)

  !! Compute the inverse square root matrix.
  IF (invec%process_grid%within_slice_rank .EQ. 0) THEN
     CALL ConstructMatrixDFromS(local_sparse, local_dense)
     CALL EigenDecomposition(local_dense, local_v, local_w)
     CALL TransposeMatrix(local_v, local_vt)
     IF (invec%is_complex) THEN
        CALL ConjugateMatrix(local_vt)
     END IF
     DO II = 1, num_vecs
        temp_val = 1.0_NTREAL/SQRT(local_w%data(II,1))
        local_vt%data(II,:) = temp_val * local_vt%data(II,:)
     END DO
     CALL MultiplyMatrix(local_v, local_vt, local_dense)
     CALL ConstructMatrixSFromD(local_dense, local_sparse)
     CALL MatrixToTripletList(local_sparse, triplet_list)
  ELSE
     CALL ConstructTripletList(triplet_list)
  END IF

  !! Rotate
  CALL ConstructEmptyMatrix(rotation, invec)
  CALL FillMatrixFromTripletList(rotation, triplet_list, &
       & preduplicated_in=.TRUE.)
  CALL MatrixMultiply(invec, rotation, outvec, memory_pool_in=pool)

  !! Cleanup
  CALL DestructMatrix(local_sparse)
  CALL DestructMatrix(local_dense)
  CALL DestructMatrix(local_v)
  CALL DestructMatrix(local_s)
  CALL DestructMatrix(local_vt)
  CALL DestructTripletList(triplet_list)
  CALL DestructMatrix(invecT)
  CALL DestructMatrix(S)
  CALL DestructMatrix(rotation)
  CALL DestructMatrix(local_w)
