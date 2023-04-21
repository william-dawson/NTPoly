  !! Local Variables
  INTEGER :: err
  INTEGER :: II
  INTEGER :: inner_len_j

  !! Local Dot
  !$omp parallel private(inner_len_j)
  !$omp do
  DO II = 1, SIZE(num_values_j)
     inner_len_j = num_values_j(II)
     out_values(II) = DotSparseVectors(indices_i(:num_values_i), &
          & values_i(:num_values_i), indices_j(:inner_len_j, II), &
          & values_j(:inner_len_j, II))
  END DO
  !$omp end do
  !$omp end parallel

  !! Reduce Over Processes
  CALL MPI_Allreduce(MPI_IN_PLACE, out_values, SIZE(num_values_j), &
       & MPINTREAL, MPI_SUM, comm, err)
