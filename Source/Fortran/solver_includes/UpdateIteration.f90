  !! Local variables
  TYPE(Matrix_ps) :: vt
  TYPE(Matrix_ps) :: Av
  TYPE(Matrix_ps) :: vAv
  TYPE(Matrix_ps) :: rotation

  !! Compute the overlap matrix.
  CALL MatrixMultiply(mat, invec, Av, memory_pool_in=pool)
  CALL TransposeMatrix(invec, vt)
  IF (mat%is_complex) THEN
     CALL ConjugateMatrix(vt)
  END IF
  CALL MatrixMultiply(vt, Av, vAv, memory_pool_in=pool)

  !! Slice out the reduced matrix.
  CALL ResizeMatrix(vAv, num_vecs)
  CALL GatherMatrixToProcess(vAv, local_sparse, 0)

  !! Compute the update
  IF (invec%process_grid%within_slice_rank .EQ. 0) THEN
     CALL ConstructMatrixDFromS(local_sparse, local_dense)

     !! Compute the eigendecomposition
     CALL EigenDecomposition(local_dense, local_v, ritz_values)

     !! Rotate
     CALL ConstructMatrixSFromD(local_v, local_sparse)
     CALL MatrixToTripletList(local_sparse, triplet_list)
  ELSE
     CALL ConstructTripletList(triplet_list)
     CALL ConstructEmptyMatrix(ritz_values, num_vecs, 1)
     ritz_values%data = 0
  END IF
  CALL ConstructEmptyMatrix(rotation, invec)
  CALL FillMatrixFromTripletList(rotation, triplet_list, &
       & preduplicated_in=.TRUE.)
  CALL MatrixMultiply(invec, rotation, outvec)

  !! Cleanup
  CALL DestructMatrix(local_sparse)
  CALL DestructMatrix(local_dense)
  CALL DestructMatrix(local_v)
  CALL DestructTripletList(triplet_list)
  CALL DestructMatrix(vt)
  CALL DestructMatrix(Av)
  CALL DestructMatrix(vAv)
  CALL DestructMatrix(rotation)
