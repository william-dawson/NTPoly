! The same NTPoly TRS2 driver built either without accelerator flags or
! against the OpenACC library. No cuSPARSE/CUDA wrapper is linked here.
PROGRAM TRS2Reference
  USE MPI
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  USE DataTypesModule, ONLY: NTREAL
  USE ProcessGridModule
  USE PSMatrixModule
  USE PSMatrixAlgebraModule, ONLY: MatrixTrace, MatrixGrandSum
  USE DensityMatrixSolversModule, ONLY: TRS2
  USE EigenBoundsModule, ONLY: GershgorinBounds
  USE SolverParametersModule
  USE PermutationModule
  USE LoadBalancerModule, ONLY: PermuteMatrix
  USE LoggingModule, ONLY: ActivateLogger, DeactivateLogger
  USE TripletModule, ONLY: Triplet_r
  USE TripletListModule, ONLY: TripletList_r, ConstructTripletList, AppendToTripletList, DestructTripletList
#ifdef REFERENCE_OPENACC
  USE openacc
#endif
  IMPLICIT NONE
  TYPE(Matrix_ps) :: input, h, identity, density
  TYPE(Permutation_t) :: permutation
  TYPE(SolverParameters_t) :: params
  CHARACTER(1024) :: filename, arg, output_density
  CHARACTER(16) :: basis_mode
  TYPE(TripletList_r) :: basis_entries
  TYPE(Triplet_r) :: basis_entry
  INTEGER :: rank, nranks, ierr, provided, i, seed_size, reps, prows, pcols
  INTEGER :: max_iter, shared, local_rank, local_size
  INTEGER, ALLOCATABLE :: seed(:)
  REAL(NTREAL) :: nelec, threshold, tolerance, energy, trace, grandsum
  REAL(NTREAL) :: emin, emax
  REAL(8) :: t0, elapsed, slowest
  CALL MPI_Init_thread(MPI_THREAD_SERIALIZED, provided, ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)
  ! Positional arguments keep the benchmark configuration explicit.
  ! input prows pcols nelec threshold tolerance max_iter repetitions
  IF (COMMAND_ARGUMENT_COUNT() < 8 .OR. COMMAND_ARGUMENT_COUNT() > 10) THEN
    IF (rank==0) PRINT *, 'input prows pcols nelec threshold tolerance max_iter repetitions [output|-] [identity|banded]'
    CALL MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  END IF
  CALL GET_COMMAND_ARGUMENT(1,filename)
  CALL GET_COMMAND_ARGUMENT(2,arg); READ(arg,*) prows
  CALL GET_COMMAND_ARGUMENT(3,arg); READ(arg,*) pcols
  CALL GET_COMMAND_ARGUMENT(4,arg); READ(arg,*) nelec
  CALL GET_COMMAND_ARGUMENT(5,arg); READ(arg,*) threshold
  CALL GET_COMMAND_ARGUMENT(6,arg); READ(arg,*) tolerance
  CALL GET_COMMAND_ARGUMENT(7,arg); READ(arg,*) max_iter
  CALL GET_COMMAND_ARGUMENT(8,arg); READ(arg,*) reps
  output_density = ''
  IF (COMMAND_ARGUMENT_COUNT()>=9) CALL GET_COMMAND_ARGUMENT(9,output_density)
  IF(TRIM(output_density)=='-') output_density=''
  basis_mode='identity'
  IF(COMMAND_ARGUMENT_COUNT()==10) CALL GET_COMMAND_ARGUMENT(10,basis_mode)
  IF(TRIM(basis_mode)/='identity' .AND. TRIM(basis_mode)/='banded') CALL MPI_Abort(MPI_COMM_WORLD,4,ierr)
  IF (prows*pcols/=nranks .OR. reps<1) CALL MPI_Abort(MPI_COMM_WORLD,2,ierr)
#ifdef REFERENCE_OPENACC
  CALL MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, shared, ierr)
  CALL MPI_Comm_rank(shared,local_rank,ierr)
  CALL MPI_Comm_size(shared,local_size,ierr)
  CALL ACC_SET_DEVICE_NUM(local_rank, ACC_DEVICE_NVIDIA)
  CALL MPI_Comm_free(shared,ierr)
  IF (rank==0) PRINT *, 'backend: OpenACC'
#else
#ifdef NTPOLY_ENABLE_GPU_TRS2
  IF (rank==0) PRINT *, 'backend: NTPoly library TRS2 (NTPOLY_GPU_TRS2 selects GPU iteration)'
#else
  IF (rank==0) PRINT *, 'backend: CPU (OpenACC directives disabled at compile time)'
#endif
#endif
  CALL ConstructProcessGrid(MPI_COMM_WORLD,prows,pcols,1)
  CALL ConstructMatrixFromBinary(input,TRIM(filename))
  CALL RANDOM_SEED(SIZE=seed_size)
  ALLOCATE(seed(seed_size))
  DO i=1,seed_size
    seed(i)=424242+i
  END DO
  CALL RANDOM_SEED(PUT=seed)
  CALL ConstructRandomPermutation(permutation,input%actual_matrix_dimension)
  CALL PermuteMatrix(input,h,permutation)
  CALL DestructPermutation(permutation)
  CALL DestructMatrix(input)
  CALL GershgorinBounds(h,emin,emax)
  IF(rank==0) THEN
    WRITE(*,*) 'actual/logical dimensions:',h%actual_matrix_dimension,h%logical_matrix_dimension
    WRITE(*,*) 'local block grid:',h%process_grid%number_of_blocks_rows,h%process_grid%number_of_blocks_columns
    WRITE(*,*) 'Gershgorin bounds:',emin,emax
  END IF
  CALL ConstructEmptyMatrix(identity,h)
  IF(TRIM(basis_mode)=='identity') THEN
    CALL FillMatrixIdentity(identity)
  ELSE
    ! A positive-definite tridiagonal inverse-overlap factor: I + 0.01*T.
    ! Constructed in the working input ordering, outside timed TRS2 calls.
    CALL ConstructTripletList(basis_entries)
    IF(rank==0) THEN
      DO i=1,h%actual_matrix_dimension
        basis_entry%index_row=i; basis_entry%index_column=i; basis_entry%point_value=1.0_NTREAL
        CALL AppendToTripletList(basis_entries,basis_entry)
        IF(i<h%actual_matrix_dimension) THEN
          basis_entry%index_row=i; basis_entry%index_column=i+1; basis_entry%point_value=0.01_NTREAL
          CALL AppendToTripletList(basis_entries,basis_entry)
          basis_entry%index_row=i+1; basis_entry%index_column=i
          CALL AppendToTripletList(basis_entries,basis_entry)
        END IF
      END DO
    END IF
    CALL FillMatrixFromTripletList(identity,basis_entries)
    CALL DestructTripletList(basis_entries)
  END IF
  IF(rank==0) PRINT *, 'reference inverse-overlap factor:',TRIM(basis_mode)
  CALL ConstructSolverParameters(params,converge_diff_in=tolerance, &
       threshold_in=threshold,max_iterations_in=max_iter,be_verbose_in=.TRUE.)
  DO i=1,reps
    CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
    IF(rank==0) THEN
      PRINT *, 'reference repetition:',i
      CALL ActivateLogger()
    END IF
    t0=MPI_Wtime()
    CALL TRS2(h,identity,nelec,density,energy_value_out=energy,solver_parameters_in=params)
    elapsed=MPI_Wtime()-t0
    IF(rank==0) CALL DeactivateLogger()
    CALL MPI_Allreduce(elapsed,slowest,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
    CALL MatrixTrace(density,trace)
    CALL MatrixGrandSum(density,grandsum)
    IF (.NOT. IEEE_IS_FINITE(energy) .OR. .NOT. IEEE_IS_FINITE(trace) .OR. &
         .NOT. IEEE_IS_FINITE(grandsum)) THEN
      IF(rank==0) PRINT *, 'ERROR: reference produced a nonfinite result'
      CALL MPI_Abort(MPI_COMM_WORLD,3,ierr)
    END IF
    IF(rank==0) THEN
      WRITE(*,'(A,I0,4(1X,ES24.16))') 'reference_result ',i,slowest,energy,trace,grandsum
      FLUSH(6)
    END IF
    IF (i==reps .AND. LEN_TRIM(output_density)>0) &
         CALL WriteMatrixToBinary(density,TRIM(output_density))
    CALL DestructMatrix(density)
  END DO
  CALL DestructSolverParameters(params)
  CALL DestructMatrix(identity)
  CALL DestructMatrix(h)
  CALL DestructProcessGrid
  CALL MPI_Finalize(ierr)
END PROGRAM
