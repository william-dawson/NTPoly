!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Solving Quantum Chemistry Systems using Purification.
MODULE DensityMatrixSolversModule
  USE DataTypesModule
  USE DistributedMatrixMemoryPoolModule
  USE DistributedSparseMatrixModule
  USE EigenBoundsModule
  USE IterativeSolversModule
  USE LoadBalancerModule
  USE LoggingModule
  USE ProcessGridModule
  USE TimerModule
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Solvers
  PUBLIC :: TRS2
  PUBLIC :: TRS4
  PUBLIC :: HPCP
  PUBLIC :: HPCPPlus
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the TRS2 method.
  !! Based on the TRS2 algorithm presented in \cite niklasson2002
  !! @param[in] Hamiltonian the matrix to compute the corresponding density from
  !! @param[in] InverseSquareRoot of the overlap matrix.
  !! @param[in] nel the number of electrons.
  !! @param[out] Density the density matrix computed by this routine.
  !! @param[out] chemical_potential_out the chemical potential calculated.
  !! @param[in] solver_parameters_in parameters for the solver
  SUBROUTINE TRS2(Hamiltonian, InverseSquareRoot, nel, Density, &
       & chemical_potential_out, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in)  :: Hamiltonian, InverseSquareRoot
    INTEGER, INTENT(in) :: nel
    TYPE(DistributedSparseMatrix), INTENT(inout) :: Density
    REAL(NTREAL), INTENT(out), OPTIONAL :: chemical_potential_out
    TYPE(IterativeSolverParameters), INTENT(in), OPTIONAL :: solver_parameters_in
    REAL(NTREAL), PARAMETER :: TWO = 2.0
    REAL(NTREAL), PARAMETER :: NEGATIVE_ONE = -1.0
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters) :: solver_parameters
    !! Local Matrices
    TYPE(DistributedSparseMatrix) :: WorkingHamiltonian
    TYPE(DistributedSparseMatrix) :: Identity
    TYPE(DistributedSparseMatrix) :: X_k, X_k2, TempMat
    !! Local Variables
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: sigma_array
    REAL(NTREAL) :: trace_value
    REAL(NTREAL) :: last_trace_value
    REAL(NTREAL) :: norm_value
    !! For computing the chemical potential
    REAL(NTREAL) :: zero_value, midpoint, interval_a, interval_b
    !! Temporary Variables
    TYPE(DistributedMatrixMemoryPool_t) :: pool1
    INTEGER :: outer_counter, inner_counter
    INTEGER :: min_size, max_size
    REAL(NTREAL) :: sparsity
    INTEGER :: total_iterations

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Density Matrix Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="TRS2")
       CALL WriteCitation("niklasson2002expansion")
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    ALLOCATE(sigma_array(solver_parameters%max_iterations))

    !! Construct All The Necessary Matrices
    CALL ConstructEmpty(Density, Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmpty(WorkingHamiltonian, Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmpty(X_k, Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmpty(X_k2, Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmpty(TempMat, Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmpty(Identity, Hamiltonian%actual_matrix_dimension)
    CALL FillDistributedIdentity(Identity)

    !! Compute the working hamiltonian.
    CALL DistributedGemm(InverseSquareRoot,Hamiltonian,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
    CALL DistributedGemm(TempMat,InverseSquareRoot,WorkingHamiltonian, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(WorkingHamiltonian, WorkingHamiltonian, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
    END IF

    !! Compute the lambda scaling value.
    CALL GershgorinBounds(WorkingHamiltonian,e_min,e_max)

    !! Initialize
    CALL CopyDistributedSparseMatrix(WorkingHamiltonian,X_k)
    CALL ScaleDistributedSparseMatrix(X_k,REAL(-1.0,NTREAL))
    CALL IncrementDistributedSparseMatrix(Identity,X_k,alpha_in=e_max)
    CALL ScaleDistributedSparseMatrix(X_k,REAL(1.0,NTREAL)/(e_max-e_min))
    trace_value = 0.0d+0
    last_trace_value = 0.0d+0

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0d+0
    DO outer_counter = 1,solver_parameters%max_iterations
       !do outer_counter = 1,1
       IF (solver_parameters%be_verbose .AND. outer_counter .GT. 1) THEN
          CALL WriteListElement(key="Round", int_value_in=outer_counter-1)
          CALL EnterSubLog
          CALL WriteListElement(key="Convergence", float_value_in=norm_value)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF

       !! Compute Sigma, as well as convergence check
       !last_trace_value = trace_value
       trace_value = Trace(X_k)
       !norm_value = abs(trace_value - last_trace_value)
       IF (nel*0.5 - trace_value .LT. REAL(0.0,NTREAL)) THEN
          sigma_array(outer_counter) = REAL(-1.0,NTREAL)
       ELSE
          sigma_array(outer_counter) = REAL(1.0,NTREAL)
       END IF

       !! Compute X_k2
       CALL DistributedGemm(X_k,X_k,X_k2, &
            & threshold_in=solver_parameters%threshold, &
            & memory_pool_in=pool1)

       !! Subtract from X_k
       CALL CopyDistributedSparseMatrix(X_k,TempMat)
       CALL IncrementDistributedSparseMatrix(X_k2,TempMat,REAL(-1.0,NTREAL))

       !! Add To X_k
       CALL CopyDistributedSparseMatrix(X_k,X_k2)
       CALL IncrementDistributedSparseMatrix(TempMat,X_k, &
            & sigma_array(outer_counter))
       CALL IncrementDistributedSparseMatrix(X_k,X_k2,alpha_in=NEGATIVE_ONE)
       norm_value = ABS(DistributedSparseNorm(X_k2))
    END DO
    total_iterations = outer_counter-1
    IF (solver_parameters%be_verbose) THEN
       CALL GetLoadBalance(X_k,min_size,max_size)
       sparsity = &
            & REAL(GetSize(X_k),NTREAL)/X_k%actual_matrix_dimension**2
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations",int_value_in=total_iterations)
       CALL WriteHeader("Load_Balance")
       CALL EnterSubLog
       CALL WriteListElement(key="min_size", int_value_in=min_size)
       CALL WriteListElement(key="max_size", int_value_in=max_size)
       CALL ExitSubLog
       CALL WriteElement(key="Sparsity", float_value_in=sparsity)
    END IF

    !! Undo Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(X_k, X_k, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
    END IF

    !! Compute the density matrix in the non-orthogonalized basis
    CALL DistributedGemm(InverseSquareRoot,X_k,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
    CALL DistributedGemm(TempMat,InverseSquareRoot,Density, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)

    !! Cleanup
    CALL DestructDistributedSparseMatrix(WorkingHamiltonian)
    CALL DestructDistributedSparseMatrix(Identity)
    CALL DestructDistributedSparseMatrix(TempMat)
    CALL DestructDistributedSparseMatrix(X_k)
    CALL DestructDistributedSparseMatrix(X_k2)
    CALL DestructDistributedMatrixMemoryPool(pool1)

    !! Compute The Chemical Potential
    IF (PRESENT(chemical_potential_out)) THEN
       interval_a = 0.0d+0
       interval_b = 1.0d+0
       midpoint = 0.0d+0
       midpoints: DO outer_counter = 1,solver_parameters%max_iterations
          midpoint = (interval_b - interval_a)/TWO + interval_a
          zero_value = midpoint
          !! Compute polynomial function at the guess point.
          polynomial:DO inner_counter=1,total_iterations
             IF (sigma_array(inner_counter) .LT. 0) THEN
                zero_value = zero_value*zero_value
             ELSE
                zero_value = 2*zero_value - zero_value*zero_value
             END IF
          END DO polynomial
          !! Change bracketing.
          IF (zero_value .LT. 0.5) THEN
             interval_a = midpoint
          ELSE
             interval_b = midpoint
          END IF
          !! Check convergence.
          IF ( ABS(zero_value-0.5) .LT. solver_parameters%converge_diff ) THEN
             EXIT
          END IF
       END DO midpoints
       !! Undo scaling.
       chemical_potential_out = e_max + (e_min - e_max)*midpoint
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    !! Cleanup
    DEALLOCATE(sigma_array)
  END SUBROUTINE TRS2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the TRS4 method.
  !! Based on the TRS4 algorithm presented in \cite niklasson2002
  !! @param[in] Hamiltonian the matrix to compute the corresponding density from
  !! @param[in] InverseSquareRoot of the overlap matrix.
  !! @param[in] nel the number of electrons.
  !! @param[out] Density the density matrix computed by this routine.
  !! @param[out] chemical_potential_out the chemical potential calculated.
  !! @param[in] solver_parameters_in parameters for the solver
  SUBROUTINE TRS4(Hamiltonian, InverseSquareRoot, nel, Density, &
       & chemical_potential_out, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in)  :: Hamiltonian, InverseSquareRoot
    INTEGER, INTENT(in) :: nel
    TYPE(DistributedSparseMatrix), INTENT(inout) :: Density
    REAL(NTREAL), INTENT(out), OPTIONAL :: chemical_potential_out
    TYPE(IterativeSolverParameters), INTENT(in), OPTIONAL :: solver_parameters_in
    REAL(NTREAL), PARAMETER :: sigma_min = 0.0
    REAL(NTREAL), PARAMETER :: sigma_max = 6.0
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters) :: solver_parameters
    !! Local Matrices
    TYPE(DistributedSparseMatrix) :: WorkingHamiltonian
    TYPE(DistributedSparseMatrix) :: Identity
    TYPE(DistributedSparseMatrix) :: X_k, X_k2, Fx_right, GX_right, TempMat
    !! Local Variables
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: sigma_array
    REAL(NTREAL) :: trace_value
    REAL(NTREAL) :: last_trace_value
    REAL(NTREAL) :: norm_value
    !! For computing the chemical potential
    REAL(NTREAL) :: zero_value, midpoint, interval_a, interval_b
    REAL(NTREAL) :: tempfx,tempgx
    !! Temporary Variables
    TYPE(DistributedMatrixMemoryPool_t) :: pool1
    INTEGER :: outer_counter, inner_counter
    INTEGER :: total_iterations
    INTEGER :: min_size, max_size
    REAL(NTREAL) :: sparsity
    REAL(NTREAL) :: trace_fx, trace_gx

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL EnterSubLog
       CALL WriteHeader("Density Matrix Solver")
       CALL WriteElement(key="Method", text_value_in="TRS4")
       CALL WriteCitation("niklasson2002expansion")
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    ALLOCATE(sigma_array(solver_parameters%max_iterations))

    !! Construct All The Necessary Matrices
    CALL ConstructEmpty(Density, Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmpty(WorkingHamiltonian, Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmpty(X_k, Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmpty(X_k2, Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmpty(TempMat, Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmpty(Fx_right, Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmpty(Gx_right, Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmpty(Identity, Hamiltonian%actual_matrix_dimension)
    CALL FillDistributedIdentity(Identity)

    !! Compute the working hamiltonian.
    CALL DistributedGemm(InverseSquareRoot,Hamiltonian,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
    CALL DistributedGemm(TempMat,InverseSquareRoot,WorkingHamiltonian, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(WorkingHamiltonian, WorkingHamiltonian, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
    END IF

    !! Compute the lambda scaling value.
    CALL GershgorinBounds(WorkingHamiltonian,e_min,e_max)

    !! Initialize
    CALL CopyDistributedSparseMatrix(WorkingHamiltonian,X_k)
    CALL ScaleDistributedSparseMatrix(X_k,REAL(-1.0,NTREAL))
    CALL IncrementDistributedSparseMatrix(Identity,X_k,alpha_in=e_max)
    CALL ScaleDistributedSparseMatrix(X_k,REAL(1.0,NTREAL)/(e_max-e_min))
    trace_value = 0.0d+0
    last_trace_value = 0.0d+0

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0d+0
    DO outer_counter = 1,solver_parameters%max_iterations
       IF (solver_parameters%be_verbose .AND. outer_counter .GT. 1) THEN
          CALL WriteListElement(key="Round", int_value_in=outer_counter-1)
          CALL EnterSubLog
          CALL WriteListElement(key="Convergence", float_value_in=norm_value)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF

       !! Compute X_k2
       CALL DistributedGemm(X_k, X_k, X_k2, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
       !! Compute Fx_right
       CALL CopyDistributedSparseMatrix(X_k2,Fx_right)
       CALL ScaleDistributedSparseMatrix(Fx_right,REAL(-3.0,NTREAL))
       CALL IncrementDistributedSparseMatrix(X_k,Fx_right, &
            & alpha_in=REAL(4.0,NTREAL))
       !! Compute Gx_right
       CALL CopyDistributedSparseMatrix(Identity,Gx_right)
       CALL IncrementDistributedSparseMatrix(X_k,Gx_right, &
            & alpha_in=REAL(-2.0,NTREAL))
       CALL IncrementDistributedSparseMatrix(X_k2,Gx_right)

       !! Compute Traces
       trace_fx = DotDistributedSparseMatrix(X_k2,Fx_right)
       trace_gx = DotDistributedSparseMatrix(X_k2,Gx_right)

       !! Compute Sigma
       sigma_array(outer_counter) = (nel*0.5-trace_fx)/trace_gx

       !! Update The Matrix
       IF (sigma_array(outer_counter) .GT. sigma_max) THEN
          CALL CopyDistributedSparseMatrix(X_k, TempMat)
          CALL ScaleDistributedSparseMatrix(TempMat, REAL(2.0,NTREAL))
          CALL IncrementDistributedSparseMatrix(X_k2, TempMat, &
               & alpha_in=REAL(-1.0,NTREAL))
       ELSE IF (sigma_array(outer_counter) .LT. sigma_min) THEN
          CALL CopyDistributedSparseMatrix(X_k2, TempMat)
       ELSE
          CALL ScaleDistributedSparseMatrix(Gx_right,sigma_array(outer_counter))
          CALL IncrementDistributedSparseMatrix(Fx_right,Gx_right)
          CALL DistributedGemm(X_k2, Gx_right, TempMat, &
               & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
       END IF

       !! Check Convergence
       CALL IncrementDistributedSparseMatrix(TempMat,X_k, &
            & alpha_in=REAL(-1.0,NTREAL))
       norm_value = ABS(DistributedSparseNorm(X_k))
       CALL CopyDistributedSparseMatrix(TempMat,X_k)

       trace_value = Trace(X_k)
    END DO
    total_iterations = outer_counter-1
    IF (solver_parameters%be_verbose) THEN
       CALL GetLoadBalance(X_k,min_size,max_size)
       sparsity = &
            & REAL(GetSize(X_k),NTREAL)/X_k%actual_matrix_dimension**2
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations",int_value_in=total_iterations)
       CALL WriteHeader("Load_Balance")
       CALL EnterSubLog
       CALL WriteListElement(key="min_size", int_value_in=min_size)
       CALL WriteListElement(key="max_size", int_value_in=max_size)
       CALL ExitSubLog
       CALL WriteElement(key="Sparsity", float_value_in=sparsity)
    END IF

    !! Undo Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(X_k, X_k, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
    END IF

    !! Compute the density matrix in the non-orthogonalized basis
    CALL DistributedGemm(InverseSquareRoot,X_k,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
    CALL DistributedGemm(TempMat,InverseSquareRoot,Density, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)

    !! Cleanup
    CALL DestructDistributedSparseMatrix(WorkingHamiltonian)
    CALL DestructDistributedSparseMatrix(Identity)
    CALL DestructDistributedSparseMatrix(TempMat)
    CALL DestructDistributedSparseMatrix(X_k)
    CALL DestructDistributedSparseMatrix(X_k2)
    CALL DestructDistributedSparseMatrix(Fx_right)
    CALL DestructDistributedSparseMatrix(Gx_right)
    CALL DestructDistributedMatrixMemoryPool(pool1)

    !! Compute The Chemical Potential
    IF (PRESENT(chemical_potential_out)) THEN
       interval_a = 0.0d+0
       interval_b = 1.0d+0
       midpoint = 0.0d+0
       midpoints: DO outer_counter = 1,solver_parameters%max_iterations
          midpoint = (interval_b - interval_a)/2.0 + interval_a
          zero_value = midpoint
          !! Compute polynomial function at the guess point.
          polynomial:DO inner_counter=1,total_iterations
             IF (sigma_array(inner_counter) .GT. sigma_max) THEN
                zero_value = 2*zero_value - zero_value*zero_value
             ELSE IF (sigma_array(inner_counter) .LT. sigma_min) THEN
                zero_value = zero_value*zero_value
             ELSE
                tempfx = (zero_value*zero_value) * &
                     (4*zero_value - 3*zero_value*zero_value)
                tempgx = (zero_value*zero_value) * (1-zero_value)*(1-zero_value)
                zero_value = tempfx + sigma_array(inner_counter)*tempgx
             END IF
          END DO polynomial
          !! Change bracketing.
          IF (zero_value .LT. 0.5) THEN
             interval_a = midpoint
          ELSE
             interval_b = midpoint
          END IF
          !! Check convergence.
          IF ( ABS(zero_value-0.5) .LT. solver_parameters%converge_diff ) THEN
             EXIT
          END IF
       END DO midpoints
       !! Undo scaling.
       chemical_potential_out = e_max + (e_min - e_max)*midpoint
    END IF

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    DEALLOCATE(sigma_array)
  END SUBROUTINE TRS4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the HPCP method.
  !! Based on the algorithm presented in \cite truflandier2016communication
  !! @param[in] Hamiltonian the matrix to compute the corresponding density from.
  !! @param[in] InverseSquareRoot of the overlap matrix.
  !! @param[in] nel the number of electrons.
  !! @param[out] Density the density matrix computed by this routine.
  !! @param[out] chemical_potential_out the chemical potential calculated.
  !! @param[in] solver_parameters_in parameters for the solver
  SUBROUTINE HPCP(Hamiltonian, InverseSquareRoot, nel, Density, &
       & chemical_potential_out, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in)  :: Hamiltonian, InverseSquareRoot
    INTEGER, INTENT(in) :: nel
    TYPE(DistributedSparseMatrix), INTENT(inout) :: Density
    REAL(NTREAL), INTENT(out), OPTIONAL :: chemical_potential_out
    TYPE(IterativeSolverParameters), INTENT(in), OPTIONAL :: solver_parameters_in
    REAL(NTREAL), PARAMETER :: TWO = 2.0
    REAL(NTREAL), PARAMETER :: NEGATIVE_ONE = -1.0
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters) :: solver_parameters
    !! Local Matrices
    TYPE(DistributedSparseMatrix) :: WorkingHamiltonian
    TYPE(DistributedSparseMatrix) :: TempMat
    TYPE(DistributedSparseMatrix) :: Identity
    TYPE(DistributedSparseMatrix) :: D1, DH, DDH, D2DH
    !! Local Variables
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL) :: beta_1, beta_2
    REAL(NTREAL) :: beta, beta_bar
    REAL(NTREAL) :: sigma, sigma_bar
    REAL(NTREAL) :: mu
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: sigma_array
    REAL(NTREAL) :: trace_value
    REAL(NTREAL) :: norm_value, norm_value2
    !! For computing the chemical potential
    REAL(NTREAL) :: zero_value, midpoint, interval_a, interval_b
    !! Temporary Variables
    TYPE(DistributedMatrixMemoryPool_t) :: pool1
    INTEGER :: outer_counter, inner_counter
    INTEGER :: total_iterations
    INTEGER :: min_size, max_size
    REAL(NTREAL) :: sparsity
    INTEGER :: matrix_dimension

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Density Matrix Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="HPCP")
       CALL WriteCitation("truflandier2016communication")
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    ALLOCATE(sigma_array(solver_parameters%max_iterations))

    matrix_dimension = Hamiltonian%actual_matrix_dimension

    !! Construct All The Necessary Matrices
    CALL ConstructEmpty(Density, matrix_dimension)
    CALL ConstructEmpty(WorkingHamiltonian, matrix_dimension)
    CALL ConstructEmpty(TempMat, matrix_dimension)
    CALL ConstructEmpty(D1, matrix_dimension)
    CALL ConstructEmpty(DH, matrix_dimension)
    CALL ConstructEmpty(DDH, matrix_dimension)
    CALL ConstructEmpty(D2DH, matrix_dimension)
    CALL ConstructEmpty(Identity, matrix_dimension)
    CALL FillDistributedIdentity(Identity)

    !! Compute the working hamiltonian.
    CALL DistributedGemm(InverseSquareRoot,Hamiltonian,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
    CALL DistributedGemm(TempMat,InverseSquareRoot,WorkingHamiltonian, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(WorkingHamiltonian, WorkingHamiltonian, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
    END IF

    !! Compute the initial matrix.
    CALL GershgorinBounds(WorkingHamiltonian,e_min,e_max)
    mu = Trace(WorkingHamiltonian)/matrix_dimension
    sigma_bar = (matrix_dimension - 0.5*nel)/matrix_dimension
    sigma = 1.0 - sigma_bar
    beta = sigma/((e_max) - mu)
    beta_bar = sigma_bar/(mu - (e_min))
    beta_1 = sigma
    beta_2 = MIN(beta,beta_bar)

    !! Initialize
    CALL CopyDistributedSparseMatrix(Identity,D1)
    CALL ScaleDistributedSparseMatrix(D1,beta_1)
    CALL CopyDistributedSparseMatrix(Identity,TempMat)
    CALL ScaleDistributedSparseMatrix(TempMat,mu)
    CALL IncrementDistributedSparseMatrix(WorkingHamiltonian, TempMat, &
         & NEGATIVE_ONE)
    CALL ScaleDistributedSparseMatrix(TempMat,beta_2)
    CALL IncrementDistributedSparseMatrix(TempMat,D1)
    trace_value = 0.0d+0

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0d+0
    norm_value2 = norm_value
    DO outer_counter = 1,solver_parameters%max_iterations
       IF (solver_parameters%be_verbose .AND. outer_counter .GT. 1) THEN
          CALL WriteListElement(key="Round", int_value_in=outer_counter-1)
          CALL EnterSubLog
          CALL WriteListElement(key="Convergence", float_value_in=norm_value)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF

       !! Compute the hole matrix DH
       CALL CopyDistributedSparseMatrix(D1,DH)
       CALL IncrementDistributedSparseMatrix(Identity,DH,alpha_in=NEGATIVE_ONE)
       CALL ScaleDistributedSparseMatrix(DH,NEGATIVE_ONE)

       !! Compute DDH, as well as convergence check
       CALL DistributedGemm(D1,DH,DDH,threshold_in=solver_parameters%threshold,&
            & memory_pool_in=pool1)
       trace_value = Trace(DDH)
       norm_value = ABS(trace_value)

       !! Compute D2DH
       CALL DistributedGemm(D1,DDH,D2DH,threshold_in=solver_parameters%threshold,&
            & memory_pool_in=pool1)

       !! Compute Sigma
       sigma_array(outer_counter) = Trace(D2DH)/trace_value

       CALL CopyDistributedSparseMatrix(D1,TempMat)

       !! Compute D1 + 2*D2DH
       CALL IncrementDistributedSparseMatrix(D2DH,D1,alpha_in=TWO)

       !! Compute D1 + 2*D2DH -2*Sigma*DDH
       CALL IncrementDistributedSparseMatrix(DDH, D1, &
            & alpha_in=NEGATIVE_ONE*TWO*sigma_array(outer_counter))

       !! Check Convergence
       CALL IncrementDistributedSparseMatrix(D1,TempMat,alpha_in=NEGATIVE_ONE)
       norm_value2 = ABS(DistributedSparseNorm(TempMat))
    END DO
    total_iterations = outer_counter-1
    IF (solver_parameters%be_verbose) THEN
       CALL GetLoadBalance(D1,min_size,max_size)
       sparsity = REAL(GetSize(D1),NTREAL)/D1%actual_matrix_dimension**2
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations",int_value_in=total_iterations)
       CALL WriteHeader("Load_Balance")
       CALL EnterSubLog
       CALL WriteListElement(key="min_size", int_value_in=min_size)
       CALL WriteListElement(key="max_size", int_value_in=max_size)
       CALL ExitSubLog
       CALL WriteElement(key="Sparsity", float_value_in=sparsity)
    END IF

    !! Undo Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(D1, D1, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
    END IF

    !! Compute the density matrix in the non-orthogonalized basis
    CALL DistributedGemm(InverseSquareRoot,D1,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
    CALL DistributedGemm(TempMat,InverseSquareRoot,Density, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)

    !! Cleanup
    CALL DestructDistributedSparseMatrix(WorkingHamiltonian)
    CALL DestructDistributedSparseMatrix(Identity)
    CALL DestructDistributedSparseMatrix(TempMat)
    CALL DestructDistributedSparseMatrix(D1)
    CALL DestructDistributedSparseMatrix(DH)
    CALL DestructDistributedSparseMatrix(DDH)
    CALL DestructDistributedSparseMatrix(D2DH)
    CALL DestructDistributedMatrixMemoryPool(pool1)

    !! Compute The Chemical Potential
    IF (PRESENT(chemical_potential_out)) THEN
       interval_a = 0.0d+0
       interval_b = 1.0d+0
       midpoint = 0.0d+0
       midpoints: DO outer_counter = 1,solver_parameters%max_iterations
          midpoint = (interval_b - interval_a)/TWO + interval_a
          zero_value = midpoint
          !! Compute polynomial function at the guess point.
          polynomial:DO inner_counter=1,total_iterations
             zero_value = zero_value + 2*((zero_value*zero_value)*(1-zero_value) - &
                  & sigma_array(inner_counter)*(zero_value*(1-zero_value)))
          END DO polynomial
          !! Change bracketing.
          IF (zero_value .LT. 0.5) THEN
             interval_a = midpoint
          ELSE
             interval_b = midpoint
          END IF
          !! Check convergence.
          IF ( ABS(zero_value-0.5) .LT. solver_parameters%converge_diff ) THEN
             EXIT
          END IF
       END DO midpoints
       !! Undo scaling.
       chemical_potential_out = mu + (beta_1 - midpoint)/beta_2
    END IF
    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF
    DEALLOCATE(sigma_array)
  END SUBROUTINE HPCP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the HPCP+ method.
  !! Based on the algorithm presented in \cite truflandier2016communication
  !! @param[in] Hamiltonian the matrix to compute the corresponding density from.
  !! @param[in] InverseSquareRoot of the overlap matrix.
  !! @param[in] nel the number of electrons.
  !! @param[out] Density the density matrix computed by this routine.
  !! @param[out] chemical_potential_out the chemical potential calculated.
  !! @param[in] solver_parameters_in parameters for the solver
  SUBROUTINE HPCPPlus(Hamiltonian, InverseSquareRoot, nel, Density, &
       & chemical_potential_out, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in)  :: Hamiltonian, InverseSquareRoot
    INTEGER, INTENT(in) :: nel
    TYPE(DistributedSparseMatrix), INTENT(inout) :: Density
    REAL(NTREAL), INTENT(out), OPTIONAL :: chemical_potential_out
    TYPE(IterativeSolverParameters), INTENT(in), OPTIONAL :: solver_parameters_in
    REAL(NTREAL), PARAMETER :: TWO = 2.0
    REAL(NTREAL), PARAMETER :: NEGATIVE_ONE = -1.0
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters) :: solver_parameters
    !! Local Matrices
    TYPE(DistributedSparseMatrix) :: WorkingHamiltonian
    TYPE(DistributedSparseMatrix) :: TempMat
    TYPE(DistributedSparseMatrix) :: Identity
    TYPE(DistributedSparseMatrix) :: D1, DH, DDH, D2DH
    !! Local Variables
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL) :: beta_1, beta_2
    REAL(NTREAL) :: beta_1_h, beta_2_h
    REAL(NTREAL) :: a, b, c, d
    REAL(NTREAL) :: mixing_interior, mixing_value
    REAL(NTREAL) :: beta, beta_bar
    REAL(NTREAL) :: sigma, sigma_bar
    REAL(NTREAL) :: mu
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: sigma_array
    REAL(NTREAL) :: trace_value
    REAL(NTREAL) :: norm_value, norm_value2
    !! For computing the chemical potential
    REAL(NTREAL) :: zero_value, midpoint, interval_a, interval_b
    !! Temporary Variables
    TYPE(DistributedMatrixMemoryPool_t) :: pool1
    INTEGER :: outer_counter, inner_counter
    INTEGER :: total_iterations
    INTEGER :: min_size, max_size
    REAL(NTREAL) :: sparsity
    INTEGER :: matrix_dimension

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Density Matrix Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="HPCP+")
       CALL WriteCitation("truflandier2016communication")
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    ALLOCATE(sigma_array(solver_parameters%max_iterations))

    matrix_dimension = Hamiltonian%actual_matrix_dimension

    !! Construct All The Necessary Matrices
    CALL ConstructEmpty(Density, matrix_dimension)
    CALL ConstructEmpty(WorkingHamiltonian, matrix_dimension)
    CALL ConstructEmpty(TempMat, matrix_dimension)
    CALL ConstructEmpty(D1, matrix_dimension)
    CALL ConstructEmpty(DH, matrix_dimension)
    CALL ConstructEmpty(DDH, matrix_dimension)
    CALL ConstructEmpty(D2DH, matrix_dimension)
    CALL ConstructEmpty(Identity, matrix_dimension)
    CALL FillDistributedIdentity(Identity)

    !! Compute the working hamiltonian.
    CALL DistributedGemm(InverseSquareRoot,Hamiltonian,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
    CALL DistributedGemm(TempMat,InverseSquareRoot,WorkingHamiltonian, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(WorkingHamiltonian, WorkingHamiltonian, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
    END IF

    !! Compute the initial matrix.
    CALL GershgorinBounds(WorkingHamiltonian,e_min,e_max)
    mu = Trace(WorkingHamiltonian)/matrix_dimension
    sigma_bar = (matrix_dimension - 0.5*nel)/matrix_dimension
    sigma = 1.0 - sigma_bar
    beta = sigma/((e_max) - mu)
    beta_bar = sigma_bar/(mu - (e_min))
    beta_1 = sigma
    beta_1_h = sigma_bar
    beta_2 = MIN(beta,beta_bar)
    beta_2_h = -1.0*MAX(beta,beta_bar)

    !! Initialize
    CALL CopyDistributedSparseMatrix(Identity,D1)
    CALL ScaleDistributedSparseMatrix(D1,beta_1)
    CALL CopyDistributedSparseMatrix(Identity,TempMat)
    CALL ScaleDistributedSparseMatrix(TempMat,mu)
    CALL IncrementDistributedSparseMatrix(WorkingHamiltonian, TempMat, &
         & NEGATIVE_ONE)
    CALL ScaleDistributedSparseMatrix(TempMat,beta_2)
    CALL IncrementDistributedSparseMatrix(TempMat,D1)

    CALL CopyDistributedSparseMatrix(Identity,DH)
    CALL ScaleDistributedSparseMatrix(DH,beta_1_h)
    CALL CopyDistributedSparseMatrix(Identity,TempMat)
    CALL ScaleDistributedSparseMatrix(TempMat,mu)
    CALL IncrementDistributedSparseMatrix(WorkingHamiltonian,TempMat, &
         & REAL(-1.0,NTREAL))
    CALL ScaleDistributedSparseMatrix(TempMat,beta_2_h)
    CALL IncrementDistributedSparseMatrix(TempMat,DH)

    CALL CopyDistributedSparseMatrix(Identity,TempMat)
    CALL IncrementDistributedSparseMatrix(DH,TempMat,REAL(-1.0,NTREAL))

    a = DotDistributedSparseMatrix(D1,D1)
    b = DotDistributedSparseMatrix(TempMat,TempMat)
    c = DotDistributedSparseMatrix(D1,TempMat)
    IF (sigma < (0.3333333333)) THEN
       d = (nel*0.5) - 0.66666666*(nel*0.5)
    ELSE
       d = (nel*0.5 - 0.666666666*(matrix_dimension - nel*0.5))
    END IF

    mixing_interior = SQRT((2*c-2*b)**2 - 4*(a+b-2*c)*(b-d))
    mixing_value = ((2*b - 2*c) + mixing_interior)/(2*(a+b - 2*c))
    IF (.NOT. mixing_value .LE. 1.0 .AND. mixing_value .GE. 0) THEN
       mixing_value = ((2*b - 2*c) - mixing_interior)/(2*(a+b - 2*c))
    ENDIF

    CALL ScaleDistributedSparseMatrix(D1,REAL(mixing_value,ntreal))
    CALL CopyDistributedSparseMatrix(Identity,TempMat)
    CALL IncrementDistributedSparseMatrix(DH,TempMat,REAL(-1.0,NTREAL))
    CALL IncrementDistributedSparseMatrix(TempMat,D1,REAL(1.0-mixing_value,ntreal))
    trace_value = 0.0d+0

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0d+0
    norm_value2 = norm_value
    DO outer_counter = 1,solver_parameters%max_iterations
       IF (solver_parameters%be_verbose .AND. outer_counter .GT. 1) THEN
          CALL WriteListElement(key="Round", int_value_in=outer_counter-1)
          CALL EnterSubLog
          CALL WriteListElement(key="Convergence", float_value_in=norm_value)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF

       !! Compute the hole matrix DH
       CALL CopyDistributedSparseMatrix(D1,DH)
       CALL IncrementDistributedSparseMatrix(Identity,DH,alpha_in=NEGATIVE_ONE)
       CALL ScaleDistributedSparseMatrix(DH,NEGATIVE_ONE)

       !! Compute DDH, as well as convergence check
       CALL DistributedGemm(D1,DH,DDH,threshold_in=solver_parameters%threshold,&
            & memory_pool_in=pool1)
       trace_value = Trace(DDH)
       norm_value = ABS(trace_value)

       !! Compute D2DH
       CALL DistributedGemm(D1,DDH,D2DH,threshold_in=solver_parameters%threshold,&
            & memory_pool_in=pool1)

       !! Compute Sigma
       sigma_array(outer_counter) = Trace(D2DH)/trace_value

       CALL CopyDistributedSparseMatrix(D1,TempMat)

       !! Compute D1 + 2*D2DH
       CALL IncrementDistributedSparseMatrix(D2DH,D1,alpha_in=TWO)

       !! Compute D1 + 2*D2DH -2*Sigma*DDH
       CALL IncrementDistributedSparseMatrix(DDH, D1, &
            & alpha_in=NEGATIVE_ONE*TWO*sigma_array(outer_counter))

       !! Check Convergence
       CALL IncrementDistributedSparseMatrix(D1,TempMat,alpha_in=NEGATIVE_ONE)
       norm_value2 = ABS(DistributedSparseNorm(TempMat))
    END DO
    total_iterations = outer_counter-1
    IF (solver_parameters%be_verbose) THEN
       CALL GetLoadBalance(D1,min_size,max_size)
       sparsity = REAL(GetSize(D1),NTREAL)/D1%actual_matrix_dimension**2
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations",int_value_in=total_iterations)
       CALL WriteHeader("Load_Balance")
       CALL EnterSubLog
       CALL WriteListElement(key="min_size", int_value_in=min_size)
       CALL WriteListElement(key="max_size", int_value_in=max_size)
       CALL ExitSubLog
       CALL WriteElement(key="Sparsity", float_value_in=sparsity)
    END IF

    !! Undo Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(D1, D1, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
    END IF

    !! Compute the density matrix in the non-orthogonalized basis
    CALL DistributedGemm(InverseSquareRoot,D1,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
    CALL DistributedGemm(TempMat,InverseSquareRoot,Density, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)

    !! Cleanup
    CALL DestructDistributedSparseMatrix(WorkingHamiltonian)
    CALL DestructDistributedSparseMatrix(Identity)
    CALL DestructDistributedSparseMatrix(TempMat)
    CALL DestructDistributedSparseMatrix(D1)
    CALL DestructDistributedSparseMatrix(DH)
    CALL DestructDistributedSparseMatrix(DDH)
    CALL DestructDistributedSparseMatrix(D2DH)
    CALL DestructDistributedMatrixMemoryPool(pool1)

    !! Compute The Chemical Potential
    IF (PRESENT(chemical_potential_out)) THEN
       interval_a = 0.0d+0
       interval_b = 1.0d+0
       midpoint = 0.0d+0
       midpoints: DO outer_counter = 1,solver_parameters%max_iterations
          midpoint = (interval_b - interval_a)/TWO + interval_a
          zero_value = midpoint
          !! Compute polynomial function at the guess point.
          polynomial:DO inner_counter=1,total_iterations
             zero_value = zero_value + 2*((zero_value*zero_value)*(1-zero_value) - &
                  & sigma_array(inner_counter)*(zero_value*(1-zero_value)))
          END DO polynomial
          !! Change bracketing.
          IF (zero_value .LT. 0.5) THEN
             interval_a = midpoint
          ELSE
             interval_b = midpoint
          END IF
          !! Check convergence.
          IF ( ABS(zero_value-0.5) .LT. solver_parameters%converge_diff ) THEN
             EXIT
          END IF
       END DO midpoints
       !! Undo scaling.
       chemical_potential_out = mu + (beta_1 - midpoint)/beta_2
    END IF
    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF
    DEALLOCATE(sigma_array)
  END SUBROUTINE HPCPPlus
END MODULE DensityMatrixSolversModule
