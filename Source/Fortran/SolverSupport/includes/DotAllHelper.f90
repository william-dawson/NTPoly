  !! Local Variables
  INTEGER :: err
  INTEGER :: counter
  INTEGER :: inner_len_j

  !! Local Dot
  !$omp parallel private(inner_len_j)
  !$omp do
  DO counter = 1, SIZE(num_values_j)
     inner_len_j = num_values_j(counter)
     out_values(counter) = DotSparseVectors(indices_i(:num_values_i), &
          & values_i(:num_values_i), indices_j(:inner_len_j, counter), &
          & values_j(:inner_len_j, counter))
  END DO
  !$omp end do
  !$omp end parallel

  !! Reduce Over Processes
  CALL MPI_Allreduce(MPI_IN_PLACE, out_values, SIZE(num_values_j), &
       & MPINTREAL, MPI_SUM, comm, err)
