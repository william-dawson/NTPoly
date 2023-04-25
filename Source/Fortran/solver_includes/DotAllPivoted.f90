  !! Local Variables
  INTEGER :: err
  INTEGER :: II
  INTEGER :: inner_len_j
  INTEGER :: local_pi_i

  !! Local Dot
  !$omp parallel private(inner_len_j, local_pi_i)
  !$omp do
  DO II = 1, num_local_pivots
     local_pi_i = pivot_vector(II)
     inner_len_j = num_values_j(local_pi_i)
     out_values(II) = DotSparseVectors(indices_i(:num_values_i), &
          & values_i(:num_values_i), indices_j(:inner_len_j, local_pi_i), &
          & values_j(:inner_len_j, local_pi_i))
  END DO
  !$omp end do
  !$omp end parallel

  !! Reduce Over Processes
  CALL MPI_Allreduce(MPI_IN_PLACE, out_values, SIZE(num_values_j), &
       & MPINTREAL, MPI_SUM, comm, err)
