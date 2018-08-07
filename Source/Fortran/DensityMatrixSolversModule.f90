!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Solving Quantum Chemistry Systems using Purification.
MODULE DensityMatrixSolversModule
  USE DataTypesModule, ONLY : NTREAL
  USE DistributedMatrixMemoryPoolModule, ONLY : DistributedMatrixMemoryPool_t, &
       & DestructDistributedMatrixMemoryPool
  USE DistributedSparseMatrixAlgebraModule, ONLY : DistributedGemm, &
       & Trace, DotDistributedSparseMatrix, IncrementDistributedSparseMatrix, &
       & ScaleDistributedSparseMatrix
  USE DistributedSparseMatrixModule, ONLY : DistributedSparseMatrix_t, &
       & ConstructEmptyDistributedSparseMatrix, CopyDistributedSparseMatrix, &
       & DestructDistributedSparseMatrix, FillDistributedIdentity, &
       & PrintMatrixInformation
  USE EigenBoundsModule, ONLY : GershgorinBounds
  USE IterativeSolversModule, ONLY : IterativeSolverParameters_t, &
       & PrintIterativeSolverParameters
  USE LoadBalancerModule, ONLY : PermuteMatrix, UndoPermuteMatrix
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteElement, &
       & WriteHeader, WriteListElement, WriteCitation
  USE ProcessGridModule
  USE TimerModule, ONLY : StartTimer, StopTimer
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Solvers
  PUBLIC :: PM
  PUBLIC :: TRS2
  PUBLIC :: TRS4
  PUBLIC :: HPCP
  ! PUBLIC :: HPCPPlus
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the PM method.
  !! Based on the PM algorithm presented in \cite palser1998
  !! @param[in] Hamiltonian the matrix to compute the corresponding density from
  !! @param[in] InverseSquareRoot of the overlap matrix.
  !! @param[in] nel the number of electrons.
  !! @param[out] Density the density matrix computed by this routine.
  !! @param[out] chemical_potential_out the chemical potential (optional).
  !! @param[in] solver_parameters_in parameters for the solver (optional)
  SUBROUTINE PM(Hamiltonian, InverseSquareRoot, nel, Density, &
       & energy_value_out, chemical_potential_out, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: Hamiltonian
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: InverseSquareRoot
    INTEGER, INTENT(IN) :: nel
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: Density
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: energy_value_out
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: chemical_potential_out
    TYPE(IterativeSolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(DistributedSparseMatrix_t) :: WorkingHamiltonian
    TYPE(DistributedSparseMatrix_t) :: Identity
    TYPE(DistributedSparseMatrix_t) :: X_k, X_k2, X_k3, TempMat
    !! Local Variables
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL) :: factor
    REAL(NTREAL) :: lambda,alpha,alpha1,alpha2
    REAL(NTREAL) :: a1,a2,a3
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: sigma_array
    REAL(NTREAL) :: trace_value
    REAL(NTREAL) :: trace_value2
    REAL(NTREAL) :: norm_value
    REAL(NTREAL) :: energy_value, energy_value2
    !! For computing the chemical potential
    REAL(NTREAL) :: zero_value, midpoint, interval_a, interval_b
    !! Temporary Variables
    TYPE(DistributedMatrixMemoryPool_t) :: pool1
    INTEGER :: outer_counter, inner_counter
    INTEGER :: total_iterations

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Density Matrix Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="PM")
       CALL WriteCitation("palser1998canonical")
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    ALLOCATE(sigma_array(solver_parameters%max_iterations))

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyDistributedSparseMatrix(Density, &
         & Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(WorkingHamiltonian, &
         & Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(X_k, &
         & Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(X_k2, &
         & Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(X_k3, &
         & Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(TempMat, &
         & Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(Identity, &
         & Hamiltonian%actual_matrix_dimension)
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

    !! Compute lambda
    trace_value = Trace(X_k)
    lambda = trace_value/X_k%actual_matrix_dimension

    !! Compute alpha
    alpha1 = nel*0.5_NTREAL/(e_max-lambda)
    alpha2 = (X_k%actual_matrix_dimension-nel*0.5_NTREAL)/(lambda-e_min)
    alpha = MIN(alpha1,alpha2)

    factor = -alpha/X_k%actual_matrix_dimension

    CALL ScaleDistributedSparseMatrix(X_k,factor)

    factor = (alpha*lambda+nel*0.5_NTREAL)/X_k%actual_matrix_dimension

    CALL IncrementDistributedSparseMatrix(Identity,X_k,alpha_in=factor)

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0_NTREAL
    energy_value = 0.0_NTREAL
    DO outer_counter = 1,solver_parameters%max_iterations
       !! Compute X_k2
       CALL DistributedGemm(X_k,X_k,X_k2, &
            & threshold_in=solver_parameters%threshold, &
            & memory_pool_in=pool1)

       !! Compute X_k3
       CALL DistributedGemm(X_k,X_k2,X_k3, &
            & threshold_in=solver_parameters%threshold, &
            & memory_pool_in=pool1)

       !! Compute X_k - X_k2
       CALL CopyDistributedSparseMatrix(X_k,TempMat)
       CALL IncrementDistributedSparseMatrix(X_k2,TempMat, &
            & alpha_in=-1.0_NTREAL, &
            & threshold_in=solver_parameters%threshold)

       !! Compute Sigma
       trace_value = Trace(TempMat)
       trace_value2 = DotDistributedSparseMatrix(TempMat,X_k)
       sigma_array(outer_counter) = trace_value2/trace_value

       IF (sigma_array(outer_counter) .GT. 0.5_NTREAL) THEN
          a1 = 0.0_NTREAL
          a2 = 1.0_NTREAL+1.0_NTREAL/sigma_array(outer_counter)
          a3 = -1.0_NTREAL/sigma_array(outer_counter)
       ELSE
          a1 = (1.0_NTREAL-2.0_NTREAL*sigma_array(outer_counter)) &
               & / (1.0_NTREAL-sigma_array(outer_counter))
          a2 = (1.0_NTREAL+sigma_array(outer_counter)) &
               & / (1.0_NTREAL-sigma_array(outer_counter))
          a3 = -1.0_NTREAL/(1.0_NTREAL-sigma_array(outer_counter))
       END IF

       !! Update X_k
       CALL ScaleDistributedSparseMatrix(X_k,a1)
       CALL IncrementDistributedSparseMatrix(X_k2,X_k, &
            & alpha_in=a2, &
            & threshold_in=solver_parameters%threshold)
       CALL IncrementDistributedSparseMatrix(X_k3,X_k, &
            & alpha_in=a3, &
            & threshold_in=solver_parameters%threshold)

       !! Energy value based convergence
       energy_value2 = energy_value
       energy_value = 2.0_NTREAL*DotDistributedSparseMatrix(X_k, &
            & WorkingHamiltonian)
       norm_value = ABS(energy_value - energy_value2)

       IF (solver_parameters%be_verbose) THEN
          CALL WriteListElement(key="Round", int_value_in=outer_counter)
          CALL EnterSubLog
          CALL WriteListElement(key="Convergence", float_value_in=norm_value)
          CALL WriteListElement("Energy_Value", float_value_in=energy_value)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    total_iterations = outer_counter-1
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations",int_value_in=outer_counter)
       CALL PrintMatrixInformation(X_k)
    END IF

    IF (PRESENT(energy_value_out)) THEN
       energy_value_out = energy_value
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
    CALL DestructDistributedSparseMatrix(X_k)
    CALL DestructDistributedSparseMatrix(X_k2)
    CALL DestructDistributedSparseMatrix(X_k3)
    CALL DestructDistributedSparseMatrix(TempMat)
    CALL DestructDistributedSparseMatrix(Identity)
    CALL DestructDistributedMatrixMemoryPool(pool1)

    !! Compute The Chemical Potential
    IF (PRESENT(chemical_potential_out)) THEN
       interval_a = 0.0_NTREAL
       interval_b = 1.0_NTREAL
       midpoint = 0.0_NTREAL
       midpoints: DO outer_counter = 1,solver_parameters%max_iterations
          midpoint = (interval_b - interval_a)/2.0_NTREAL + interval_a
          zero_value = midpoint
          !! Compute polynomial function at the guess point.
          polynomial:DO inner_counter=1,total_iterations
             IF (sigma_array(inner_counter) .GT. 0.5_NTREAL) THEN
                zero_value = ((1.0_NTREAL+ &
                     & sigma_array(inner_counter))*zero_value**2) &
                     & - (zero_value**3)
                zero_value = zero_value/sigma_array(inner_counter)
             ELSE
                zero_value = ((1.0_NTREAL - 2.0_NTREAL* &
                     & sigma_array(inner_counter))*zero_value) &
                     & + ((1.0_NTREAL+sigma_array(inner_counter))* &
                     & zero_value**2) - (zero_value**3)
                zero_value = zero_value/(1.0_NTREAL-sigma_array(inner_counter))
             END IF
          END DO polynomial
          !! Change bracketing.
          IF (zero_value .LT. 0.5_NTREAL) THEN
             interval_a = midpoint
          ELSE
             interval_b = midpoint
          END IF
          !! Check convergence.
          IF (ABS(zero_value-0.5_NTREAL) .LT. &
               & solver_parameters%converge_diff) THEN
             EXIT
          END IF
       END DO midpoints
       !! Undo scaling.
       chemical_potential_out = (nel*0.5_NTREAL-midpoint)/alpha+lambda
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    !! Cleanup
    DEALLOCATE(sigma_array)
  END SUBROUTINE PM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the TRS2 method.
  !! Based on the TRS2 algorithm presented in \cite niklasson2002
  !! @param[in] Hamiltonian the matrix to compute the corresponding density from
  !! @param[in] InverseSquareRoot of the overlap matrix.
  !! @param[in] nel the number of electrons.
  !! @param[out] Density the density matrix computed by this routine.
  !! @param[out] chemical_potential_out the chemical potential (optional).
  !! @param[in] solver_parameters_in parameters for the solver (optional)
  SUBROUTINE TRS2(Hamiltonian, InverseSquareRoot, nel, Density, &
       & energy_value_out, chemical_potential_out, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: Hamiltonian
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: InverseSquareRoot
    INTEGER, INTENT(IN) :: nel
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: Density
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: energy_value_out
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: chemical_potential_out
    TYPE(IterativeSolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(DistributedSparseMatrix_t) :: WorkingHamiltonian
    TYPE(DistributedSparseMatrix_t) :: Identity
    TYPE(DistributedSparseMatrix_t) :: X_k, X_k2, TempMat
    !! Local Variables
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: sigma_array
    REAL(NTREAL) :: trace_value
    REAL(NTREAL) :: norm_value
    REAL(NTREAL) :: energy_value, energy_value2
    !! For computing the chemical potential
    REAL(NTREAL) :: zero_value, midpoint, interval_a, interval_b
    !! Temporary Variables
    TYPE(DistributedMatrixMemoryPool_t) :: pool1
    INTEGER :: outer_counter, inner_counter
    INTEGER :: total_iterations

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters_t()
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
    CALL ConstructEmptyDistributedSparseMatrix(Density, &
         & Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(WorkingHamiltonian, &
         & Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(X_k, &
         & Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(X_k2, &
         & Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(TempMat, &
         & Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(Identity, &
         & Hamiltonian%actual_matrix_dimension)
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
    CALL ScaleDistributedSparseMatrix(X_k,-1.0_NTREAL)
    CALL IncrementDistributedSparseMatrix(Identity,X_k,alpha_in=e_max)
    CALL ScaleDistributedSparseMatrix(X_k,1.0_NTREAL/(e_max-e_min))

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0_NTREAL
    energy_value = 0.0_NTREAL
    DO outer_counter = 1,solver_parameters%max_iterations
       !! Compute Sigma
       trace_value = Trace(X_k)
       IF (nel*0.5_NTREAL - trace_value .LT. 0.0_NTREAL) THEN
          sigma_array(outer_counter) = -1.0_NTREAL
       ELSE
          sigma_array(outer_counter) = 1.0_NTREAL
       END IF

       !! Compute X_k2
       CALL DistributedGemm(X_k,X_k,X_k2, &
            & threshold_in=solver_parameters%threshold, &
            & memory_pool_in=pool1)

       !! Update X_k
       IF (sigma_array(outer_counter) .GT. 0.0_NTREAL) THEN
          CALL ScaleDistributedSparseMatrix(X_k,2.0_NTREAL)
          CALL IncrementDistributedSparseMatrix(X_k2,X_k, &
               & alpha_in=-1.0_NTREAL, &
               & threshold_in=solver_parameters%threshold)
       ELSE
          CALL CopyDistributedSparseMatrix(X_k2,X_k)
       END IF

       !! Energy value based convergence
       energy_value2 = energy_value
       energy_value = 2.0_NTREAL*DotDistributedSparseMatrix(X_k, &
            & WorkingHamiltonian)
       norm_value = ABS(energy_value - energy_value2)

       IF (solver_parameters%be_verbose) THEN
          CALL WriteListElement(key="Round", int_value_in=outer_counter)
          CALL EnterSubLog
          CALL WriteListElement(key="Convergence", float_value_in=norm_value)
          CALL WriteListElement("Energy_Value", float_value_in=energy_value)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    total_iterations = outer_counter-1
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations",int_value_in=outer_counter)
       CALL PrintMatrixInformation(X_k)
    END IF

    IF (PRESENT(energy_value_out)) THEN
       energy_value_out = energy_value
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
    CALL DestructDistributedSparseMatrix(X_k)
    CALL DestructDistributedSparseMatrix(X_k2)
    CALL DestructDistributedSparseMatrix(TempMat)
    CALL DestructDistributedSparseMatrix(Identity)
    CALL DestructDistributedMatrixMemoryPool(pool1)

    !! Compute The Chemical Potential
    IF (PRESENT(chemical_potential_out)) THEN
       interval_a = 0.0_NTREAL
       interval_b = 1.0_NTREAL
       midpoint = 0.0_NTREAL
       midpoints: DO outer_counter = 1,solver_parameters%max_iterations
          midpoint = (interval_b - interval_a)/2.0_NTREAL + interval_a
          zero_value = midpoint
          !! Compute polynomial function at the guess point.
          polynomial:DO inner_counter=1,total_iterations
             IF (sigma_array(inner_counter) .LT. 0.0_NTREAL) THEN
                zero_value = zero_value*zero_value
             ELSE
                zero_value = 2.0_NTREAL*zero_value - zero_value*zero_value
             END IF
          END DO polynomial
          !! Change bracketing.
          IF (zero_value .LT. 0.5_NTREAL) THEN
             interval_a = midpoint
          ELSE
             interval_b = midpoint
          END IF
          !! Check convergence.
          IF (ABS(zero_value-0.5_NTREAL) .LT. &
               & solver_parameters%converge_diff) THEN
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
  !! @param[out] chemical_potential_out the chemical potential (optional).
  !! @param[in] solver_parameters_in parameters for the solver (optional).
  SUBROUTINE TRS4(Hamiltonian, InverseSquareRoot, nel, Density, &
       & energy_value_out, chemical_potential_out, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN)  :: Hamiltonian
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: InverseSquareRoot
    INTEGER, INTENT(IN) :: nel
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: Density
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: energy_value_out
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: chemical_potential_out
    TYPE(IterativeSolverParameters_t), INTENT(IN), OPTIONAL :: &
         & solver_parameters_in
    REAL(NTREAL), PARAMETER :: sigma_min = 0.0_NTREAL
    REAL(NTREAL), PARAMETER :: sigma_max = 6.0_NTREAL
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(DistributedSparseMatrix_t) :: WorkingHamiltonian
    TYPE(DistributedSparseMatrix_t) :: Identity
    TYPE(DistributedSparseMatrix_t) :: X_k, X_k2, Fx_right, GX_right, TempMat
    !! Local Variables
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: sigma_array
    REAL(NTREAL) :: norm_value
    REAL(NTREAL) :: energy_value, energy_value2
    !! For computing the chemical potential
    REAL(NTREAL) :: zero_value, midpoint, interval_a, interval_b
    REAL(NTREAL) :: tempfx,tempgx
    !! Temporary Variables
    TYPE(DistributedMatrixMemoryPool_t) :: pool1
    INTEGER :: outer_counter, inner_counter
    INTEGER :: total_iterations
    REAL(NTREAL) :: trace_fx, trace_gx

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters_t()
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
    CALL ConstructEmptyDistributedSparseMatrix(Density, &
         & Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(WorkingHamiltonian, &
         & Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(X_k, &
         & Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(X_k2, &
         & Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(TempMat, &
         & Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(Fx_right, &
         & Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(Gx_right, &
         & Hamiltonian%actual_matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(Identity, &
         & Hamiltonian%actual_matrix_dimension)
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
    CALL ScaleDistributedSparseMatrix(X_k,-1.0_NTREAL)
    CALL IncrementDistributedSparseMatrix(Identity,X_k,alpha_in=e_max)
    CALL ScaleDistributedSparseMatrix(X_k,1.0_NTREAL/(e_max-e_min))

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0_NTREAL
    energy_value = 0.0_NTREAL
    DO outer_counter = 1,solver_parameters%max_iterations
       !! Compute X_k2
       CALL DistributedGemm(X_k, X_k, X_k2, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
       !! Compute Fx_right
       CALL CopyDistributedSparseMatrix(X_k2,Fx_right)
       CALL ScaleDistributedSparseMatrix(Fx_right,-3.0_NTREAL)
       CALL IncrementDistributedSparseMatrix(X_k,Fx_right, &
            & alpha_in=4.0_NTREAL)
       !! Compute Gx_right
       CALL CopyDistributedSparseMatrix(Identity,Gx_right)
       CALL IncrementDistributedSparseMatrix(X_k,Gx_right, &
            & alpha_in=-2.0_NTREAL)
       CALL IncrementDistributedSparseMatrix(X_k2,Gx_right)

       !! Compute Traces
       trace_fx = DotDistributedSparseMatrix(X_k2,Fx_right)
       trace_gx = DotDistributedSparseMatrix(X_k2,Gx_right)

       !! Avoid Overflow
       IF (ABS(trace_gx) .LT. 1.0e-14_NTREAL) THEN
          EXIT
       END IF

       !! Compute Sigma
       sigma_array(outer_counter) = (nel*0.5_NTREAL-trace_fx)/trace_gx

       !! Update The Matrix
       IF (sigma_array(outer_counter) .GT. sigma_max) THEN
          CALL CopyDistributedSparseMatrix(X_k, TempMat)
          CALL ScaleDistributedSparseMatrix(TempMat, 2.0_NTREAL)
          CALL IncrementDistributedSparseMatrix(X_k2, TempMat, &
               & alpha_in=-1.0_NTREAL)
       ELSE IF (sigma_array(outer_counter) .LT. sigma_min) THEN
          CALL CopyDistributedSparseMatrix(X_k2, TempMat)
       ELSE
          CALL ScaleDistributedSparseMatrix(Gx_right,sigma_array(outer_counter))
          CALL IncrementDistributedSparseMatrix(Fx_right,Gx_right)
          CALL DistributedGemm(X_k2, Gx_right, TempMat, &
               & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
       END IF

       CALL IncrementDistributedSparseMatrix(TempMat,X_k, &
            & alpha_in=-1.0_NTREAL)
       CALL CopyDistributedSparseMatrix(TempMat,X_k)

       !! Energy value based convergence
       energy_value2 = energy_value
       energy_value = 2.0_NTREAL*DotDistributedSparseMatrix(X_k, &
            & WorkingHamiltonian)
       norm_value = ABS(energy_value - energy_value2)

       IF (solver_parameters%be_verbose) THEN
          CALL WriteListElement(key="Round", int_value_in=outer_counter)
          CALL EnterSubLog
          CALL WriteListElement(key="Convergence", float_value_in=norm_value)
          CALL WriteListElement("Energy_Value", float_value_in=energy_value)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    total_iterations = outer_counter-1
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations",int_value_in=outer_counter)
       CALL PrintMatrixInformation(X_k)
    END IF

    IF (PRESENT(energy_value_out)) THEN
       energy_value_out = energy_value
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
    CALL DestructDistributedSparseMatrix(X_k)
    CALL DestructDistributedSparseMatrix(X_k2)
    CALL DestructDistributedSparseMatrix(Fx_right)
    CALL DestructDistributedSparseMatrix(Gx_right)
    CALL DestructDistributedSparseMatrix(TempMat)
    CALL DestructDistributedSparseMatrix(Identity)
    CALL DestructDistributedMatrixMemoryPool(pool1)

    !! Compute The Chemical Potential
    IF (PRESENT(chemical_potential_out)) THEN
       interval_a = 0.0_NTREAL
       interval_b = 1.0_NTREAL
       midpoint = 0.0_NTREAL
       midpoints: DO outer_counter = 1,solver_parameters%max_iterations
          midpoint = (interval_b - interval_a)/2.0_NTREAL + interval_a
          zero_value = midpoint
          !! Compute polynomial function at the guess point.
          polynomial:DO inner_counter=1,total_iterations
             IF (sigma_array(inner_counter) .GT. sigma_max) THEN
                zero_value = 2.0_NTREAL*zero_value - zero_value*zero_value
             ELSE IF (sigma_array(inner_counter) .LT. sigma_min) THEN
                zero_value = zero_value*zero_value
             ELSE
                tempfx = (zero_value*zero_value) &
                     & * (4.0_NTREAL*zero_value - &
                     & 3.0_NTREAL*zero_value*zero_value)
                tempgx = (zero_value*zero_value) * (1.0_NTREAL-zero_value) &
                     & * (0.0_NTREAL-zero_value)
                zero_value = tempfx + sigma_array(inner_counter)*tempgx
             END IF
          END DO polynomial
          !! Change bracketing.
          IF (zero_value .LT. 0.5_NTREAL) THEN
             interval_a = midpoint
          ELSE
             interval_b = midpoint
          END IF
          !! Check convergence.
          IF (ABS(zero_value-0.5_NTREAL) .LT. &
               & solver_parameters%converge_diff) THEN
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
  !! @param[in] Hamiltonian the matrix to compute the corresponding density from
  !! @param[in] InverseSquareRoot of the overlap matrix.
  !! @param[in] nel the number of electrons.
  !! @param[out] Density the density matrix computed by this routine.
  !! @param[out] chemical_potential_out the chemical potential (optional).
  !! @param[in] solver_parameters_in parameters for the solver (optional).
  SUBROUTINE HPCP(Hamiltonian, InverseSquareRoot, nel, Density, &
       & energy_value_out, chemical_potential_out, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: Hamiltonian
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: InverseSquareRoot
    INTEGER, INTENT(IN) :: nel
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: Density
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: energy_value_out
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: chemical_potential_out
    TYPE(IterativeSolverParameters_t), INTENT(IN), OPTIONAL :: &
         & solver_parameters_in
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(DistributedSparseMatrix_t) :: WorkingHamiltonian
    TYPE(DistributedSparseMatrix_t) :: TempMat
    TYPE(DistributedSparseMatrix_t) :: Identity
    TYPE(DistributedSparseMatrix_t) :: D1, DH, DDH, D2DH
    !! Local Variables
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL) :: beta_1, beta_2
    REAL(NTREAL) :: beta, beta_bar
    REAL(NTREAL) :: sigma, sigma_bar
    REAL(NTREAL) :: mu
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: sigma_array
    REAL(NTREAL) :: trace_value
    REAL(NTREAL) :: norm_value, norm_value2
    REAL(NTREAL) :: energy_value, energy_value2
    !! For computing the chemical potential
    REAL(NTREAL) :: zero_value, midpoint, interval_a, interval_b
    !! Temporary Variables
    TYPE(DistributedMatrixMemoryPool_t) :: pool1
    INTEGER :: outer_counter, inner_counter
    INTEGER :: total_iterations
    INTEGER :: matrix_dimension

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters_t()
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
    CALL ConstructEmptyDistributedSparseMatrix(Density, matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(WorkingHamiltonian, &
         & matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(TempMat, matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(D1, matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(DH, matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(DDH, matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(D2DH, matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(Identity, matrix_dimension)
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
    sigma_bar = (matrix_dimension - 0.5_NTREAL*nel)/matrix_dimension
    sigma = 1.0_NTREAL - sigma_bar
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
         & -1.0_NTREAL)
    CALL ScaleDistributedSparseMatrix(TempMat,beta_2)
    CALL IncrementDistributedSparseMatrix(TempMat,D1)
    trace_value = 0.0_NTREAL

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0_NTREAL
    norm_value2 = norm_value
    energy_value = 0.0_NTREAL
    DO outer_counter = 1,solver_parameters%max_iterations
       !! Compute the hole matrix DH
       CALL CopyDistributedSparseMatrix(D1,DH)
       CALL IncrementDistributedSparseMatrix(Identity,DH,alpha_in=-1.0_NTREAL)
       CALL ScaleDistributedSparseMatrix(DH,-1.0_NTREAL)

       !! Compute DDH, as well as convergence check
       CALL DistributedGemm(D1,DH,DDH,threshold_in=solver_parameters%threshold,&
            & memory_pool_in=pool1)
       trace_value = Trace(DDH)
       norm_value = ABS(trace_value)

       !! Compute D2DH
       CALL DistributedGemm(D1,DDH,D2DH, &
            & threshold_in=solver_parameters%threshold, &
            & memory_pool_in=pool1)

       !! Compute Sigma
       sigma_array(outer_counter) = Trace(D2DH)/trace_value

       CALL CopyDistributedSparseMatrix(D1,TempMat)

       !! Compute D1 + 2*D2DH
       CALL IncrementDistributedSparseMatrix(D2DH,D1,alpha_in=2.0_NTREAL)

       !! Compute D1 + 2*D2DH -2*Sigma*DDH
       CALL IncrementDistributedSparseMatrix(DDH, D1, &
            & alpha_in=-1.0_NTREAL*2.0_NTREAL*sigma_array(outer_counter))

       !! Energy value based convergence
       energy_value2 = energy_value
       energy_value = 2.0_NTREAL*DotDistributedSparseMatrix(D1, &
            & WorkingHamiltonian)
       norm_value = ABS(energy_value - energy_value2)

       IF (solver_parameters%be_verbose) THEN
          CALL WriteListElement(key="Round", int_value_in=outer_counter)
          CALL EnterSubLog
          CALL WriteListElement(key="Convergence", float_value_in=norm_value)
          CALL WriteListElement("Energy_Value", float_value_in=energy_value)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    total_iterations = outer_counter-1
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations",int_value_in=outer_counter)
       CALL PrintMatrixInformation(D1)
    END IF

    IF (PRESENT(energy_value_out)) THEN
       energy_value_out = energy_value
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
    CALL DestructDistributedSparseMatrix(TempMat)
    CALL DestructDistributedSparseMatrix(D1)
    CALL DestructDistributedSparseMatrix(DH)
    CALL DestructDistributedSparseMatrix(DDH)
    CALL DestructDistributedSparseMatrix(D2DH)
    CALL DestructDistributedSparseMatrix(Identity)
    CALL DestructDistributedMatrixMemoryPool(pool1)

    !! Compute The Chemical Potential
    IF (PRESENT(chemical_potential_out)) THEN
       interval_a = 0.0_NTREAL
       interval_b = 1.0_NTREAL
       midpoint = 0.0_NTREAL
       midpoints: DO outer_counter = 1,solver_parameters%max_iterations
          midpoint = (interval_b - interval_a)/2.0_NTREAL + interval_a
          zero_value = midpoint
          !! Compute polynomial function at the guess point.
          polynomial:DO inner_counter=1,total_iterations
             zero_value = zero_value + &
                  & 2.0_NTREAL*((zero_value**2)*(1.0_NTREAL-zero_value) &
                  & - sigma_array(inner_counter)* &
                  & zero_value*(1.0_NTREAL-zero_value))
          END DO polynomial
          !! Change bracketing.
          IF (zero_value .LT. 0.5_NTREAL) THEN
             interval_a = midpoint
          ELSE
             interval_b = midpoint
          END IF
          !! Check convergence.
          IF (ABS(zero_value-0.5_NTREAL) .LT. &
               & solver_parameters%converge_diff) THEN
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
  !! @param[out] chemical_potential_out the chemical potential (optional).
  !! @param[in] solver_parameters_in parameters for the solver (optional).
  SUBROUTINE HPCPPlus(Hamiltonian, InverseSquareRoot, nel, Density, &
       & energy_value_out, chemical_potential_out, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: Hamiltonian
    TYPE(DistributedSparseMatrix_t), INTENT(IN) :: InverseSquareRoot
    INTEGER, INTENT(IN) :: nel
    TYPE(DistributedSparseMatrix_t), INTENT(INOUT) :: Density
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: energy_value_out
    REAL(NTREAL), INTENT(OUT), OPTIONAL :: chemical_potential_out
    TYPE(IterativeSolverParameters_t), INTENT(IN), OPTIONAL :: &
         & solver_parameters_in
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(DistributedSparseMatrix_t) :: WorkingHamiltonian
    TYPE(DistributedSparseMatrix_t) :: TempMat
    TYPE(DistributedSparseMatrix_t) :: Identity
    TYPE(DistributedSparseMatrix_t) :: D1, DH, DDH, D2DH
    !! Local Variables
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL) :: beta_1, beta_2
    REAL(NTREAL) :: beta_1_h, beta_2_h
    REAL(NTREAL) :: a, b, c, d
    REAL(NTREAL) :: one_third
    REAL(NTREAL) :: mixing_interior, mixing_value
    REAL(NTREAL) :: beta, beta_bar
    REAL(NTREAL) :: sigma, sigma_bar
    REAL(NTREAL) :: mu
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: sigma_array
    REAL(NTREAL) :: trace_value
    REAL(NTREAL) :: norm_value
    REAL(NTREAL) :: energy_value, energy_value2
    !! For computing the chemical potential
    REAL(NTREAL) :: zero_value, midpoint, interval_a, interval_b
    !! Temporary Variables
    TYPE(DistributedMatrixMemoryPool_t) :: pool1
    INTEGER :: outer_counter, inner_counter
    INTEGER :: total_iterations
    INTEGER :: matrix_dimension

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters_t()
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
    CALL ConstructEmptyDistributedSparseMatrix(Density, matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(WorkingHamiltonian, &
         & matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(TempMat, matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(D1, matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(DH, matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(DDH, matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(D2DH, matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(Identity, matrix_dimension)
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
    sigma_bar = (matrix_dimension - 0.5_NTREAL*nel)/matrix_dimension
    sigma = 1.0_NTREAL - sigma_bar
    beta = sigma/((e_max) - mu)
    beta_bar = sigma_bar/(mu - (e_min))
    beta_1 = sigma
    beta_1_h = sigma_bar
    beta_2 = MIN(beta,beta_bar)
    beta_2_h = -1.0_NTREAL*MAX(beta,beta_bar)

    !! Initialize
    CALL CopyDistributedSparseMatrix(Identity,D1)
    CALL ScaleDistributedSparseMatrix(D1,beta_1)
    CALL CopyDistributedSparseMatrix(Identity,TempMat)
    CALL ScaleDistributedSparseMatrix(TempMat,mu)
    CALL IncrementDistributedSparseMatrix(WorkingHamiltonian, TempMat, &
         & -1.0_NTREAL)
    CALL ScaleDistributedSparseMatrix(TempMat,beta_2)
    CALL IncrementDistributedSparseMatrix(TempMat,D1)

    CALL CopyDistributedSparseMatrix(Identity,DH)
    CALL ScaleDistributedSparseMatrix(DH,beta_1_h)
    CALL CopyDistributedSparseMatrix(Identity,TempMat)
    CALL ScaleDistributedSparseMatrix(TempMat,mu)
    CALL IncrementDistributedSparseMatrix(WorkingHamiltonian,TempMat,-1.0_NTREAL)
    CALL ScaleDistributedSparseMatrix(TempMat,beta_2_h)
    CALL IncrementDistributedSparseMatrix(TempMat,DH)

    CALL CopyDistributedSparseMatrix(Identity,TempMat)
    CALL IncrementDistributedSparseMatrix(DH,TempMat,-1.0_NTREAL)

    a = DotDistributedSparseMatrix(D1,D1)
    b = DotDistributedSparseMatrix(TempMat,TempMat)
    c = DotDistributedSparseMatrix(D1,TempMat)

    one_third = 1.0_NTREAL/3.0_NTREAL
    IF (sigma < one_third) THEN
       d = (nel*0.5_NTREAL) - 2.0_NTREAL*one_third*(nel*0.5_NTREAL)
    ELSE
       d = (nel*0.5_NTREAL) - &
            & 2.0_NTREAL*one_third*(matrix_dimension - nel*0.5_NTREAL)
    END IF

    mixing_interior = SQRT((2.0_NTREAL*c-2.0_NTREAL*b)**2 &
         & - 4.0_NTREAL*(a+b-2.0_NTREAL*c)*(b-d))
    mixing_value = ((2.0_NTREAL*b - 2.0_NTREAL*c) + mixing_interior) &
         & / (2.0_NTREAL*(a+b - 2.0_NTREAL*c))
    IF (.NOT. mixing_value .LE. 1.0_NTREAL &
         & .AND. mixing_value .GE. 0.0_NTREAL) THEN
       mixing_value = ((2.0_NTREAL*b - 2.0_NTREAL*c) - mixing_interior) &
            & / (2.0_NTREAL*(a+b - 2.0_NTREAL*c))
    ENDIF

    CALL ScaleDistributedSparseMatrix(D1,mixing_value)
    CALL CopyDistributedSparseMatrix(Identity,TempMat)
    CALL IncrementDistributedSparseMatrix(DH,TempMat,-1.0_NTREAL)
    CALL IncrementDistributedSparseMatrix(TempMat,D1,1.0_NTREAL-mixing_value)

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0_NTREAL
    energy_value = 0.0_NTREAL
    DO outer_counter = 1,solver_parameters%max_iterations
       !! Compute the hole matrix DH
       CALL CopyDistributedSparseMatrix(D1,DH)
       CALL IncrementDistributedSparseMatrix(Identity,DH,alpha_in=-1.0_NTREAL)
       CALL ScaleDistributedSparseMatrix(DH,-1.0_NTREAL)

       !! Compute DDH, as well as convergence check
       CALL DistributedGemm(D1,DH,DDH,threshold_in=solver_parameters%threshold,&
            & memory_pool_in=pool1)
       trace_value = Trace(DDH)
       norm_value = ABS(trace_value)

       !! Compute D2DH
       CALL DistributedGemm(D1,DDH,D2DH, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)

       !! Compute Sigma
       sigma_array(outer_counter) = Trace(D2DH)/trace_value

       CALL CopyDistributedSparseMatrix(D1,TempMat)

       !! Compute D1 + 2*D2DH
       CALL IncrementDistributedSparseMatrix(D2DH,D1,alpha_in=2.0_NTREAL)

       !! Compute D1 + 2*D2DH -2*Sigma*DDH
       CALL IncrementDistributedSparseMatrix(DDH, D1, &
            & alpha_in=-1.0_NTREAL*2.0_NTREAL*sigma_array(outer_counter))

       !! Energy value based convergence
       energy_value2 = energy_value
       energy_value = 2.0_NTREAL*DotDistributedSparseMatrix(D1, &
            & WorkingHamiltonian)
       norm_value = ABS(energy_value - energy_value2)

       IF (solver_parameters%be_verbose) THEN
          CALL WriteListElement(key="Round", int_value_in=outer_counter)
          CALL EnterSubLog
          CALL WriteListElement(key="Convergence", float_value_in=norm_value)
          CALL WriteListElement("Energy_Value", float_value_in=energy_value)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    total_iterations = outer_counter-1
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations",int_value_in=outer_counter)
       CALL PrintMatrixInformation(D1)
    END IF

    IF (PRESENT(energy_value_out)) THEN
       energy_value_out = energy_value
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
    CALL DestructDistributedSparseMatrix(TempMat)
    CALL DestructDistributedSparseMatrix(D1)
    CALL DestructDistributedSparseMatrix(DH)
    CALL DestructDistributedSparseMatrix(DDH)
    CALL DestructDistributedSparseMatrix(D2DH)
    CALL DestructDistributedSparseMatrix(Identity)
    CALL DestructDistributedMatrixMemoryPool(pool1)

    !! Compute The Chemical Potential
    IF (PRESENT(chemical_potential_out)) THEN
       interval_a = 0.0_NTREAL
       interval_b = 1.0_NTREAL
       midpoint = 0.0_NTREAL
       midpoints: DO outer_counter = 1,solver_parameters%max_iterations
          midpoint = (interval_b - interval_a)/2.0_NTREAL + interval_a
          zero_value = midpoint
          !! Compute polynomial function at the guess point.
          polynomial:DO inner_counter=1,total_iterations
             zero_value = zero_value + &
                  & 2.0_NTREAL*((zero_value*zero_value)*(1.0_NTREAL-zero_value)&
                  & - sigma_array(inner_counter)* &
                  & (zero_value*(1.0_NTREAL-zero_value)))
          END DO polynomial
          !! Change bracketing.
          IF (zero_value .LT. 0.5_NTREAL) THEN
             interval_a = midpoint
          ELSE
             interval_b = midpoint
          END IF
          !! Check convergence.
          IF (ABS(zero_value-0.5_NTREAL) .LT. &
               & solver_parameters%converge_diff) THEN
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE DensityMatrixSolversModule
