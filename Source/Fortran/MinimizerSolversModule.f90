!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Solving Systems Quantum Chemistry Systems Using Minimization.
MODULE MinimizerSolversModule
  USE DataTypesModule
  USE DistributedMatrixMemoryPoolModule
  USE DistributedSparseMatrixAlgebraModule
  USE DistributedSparseMatrixModule
  USE IterativeSolversModule
  USE LoadBalancerModule
  USE LoggingModule
  USE ProcessGridModule
  USE TimerModule
  USE mpi
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Solvers
  PUBLIC :: ConjugateGradient
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the density matrix from a Hamiltonian using the CG method.
  !! Based on two papers. The first by Scuseria \cite millam1997linear developed
  !! the initial method, and then Challacombe \cite challacombe1999simplified
  !! developed a simplified scheme.
  !! @param[in] Hamiltonian the matrix to compute the corresponding density from.
  !! @param[in] InverseSquareRoot of the overlap matrix.
  !! @param[in] nel the number of electrons.
  !! @param[out] Density the density matrix computed by this routine.
  !! @param[out] chemical_potential_out the chemical potential calculated.
  !! @param[in] solver_parameters_in parameters for the solver
  !! @todo chemical potential doesn't work.
  SUBROUTINE ConjugateGradient(Hamiltonian, InverseSquareRoot, nel, Density, &
       & chemical_potential_out, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix_t), INTENT(in)  :: Hamiltonian, &
         & InverseSquareRoot
    INTEGER, INTENT(in) :: nel
    TYPE(DistributedSparseMatrix_t), INTENT(inout) :: Density
    REAL(NTREAL), INTENT(out), OPTIONAL :: chemical_potential_out
    TYPE(IterativeSolverParameters_t), INTENT(in), OPTIONAL :: &
         & solver_parameters_in
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters_t) :: solver_parameters
    !! Local Variables
    REAL(NTREAL) :: trace_value
    REAL(NTREAL) :: last_trace_value
    REAL(NTREAL) :: norm_value
    TYPE(DistributedSparseMatrix_t) :: WorkingHamiltonian
    TYPE(DistributedSparseMatrix_t) :: P_k
    TYPE(DistributedSparseMatrix_t) :: Gradient
    TYPE(DistributedSparseMatrix_t) :: G_k, G_kplusone
    TYPE(DistributedSparseMatrix_t) :: H_k
    TYPE(DistributedSparseMatrix_t) :: TempMat, TempMat2
    TYPE(DistributedSparseMatrix_t) :: Identity
    REAL(NTREAL) :: mu
    REAL(NTREAL) :: gamma
    REAL(NTREAL) :: step_size
    REAL(NTREAL) :: b, c, d
    REAL(NTREAL) :: root1, root2, root_temp
    REAL(NTREAL) :: energy_value, energy_value2
    !! Temporary Variables
    TYPE(DistributedMatrixMemoryPool_t) :: pool1
    INTEGER :: outer_counter
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
       CALL WriteElement(key="Method", text_value_in="CG")
       CALL WriteCitation("millam1997linear challacombe1999simplified")
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    !! Construct All The Necessary Matrices
    matrix_dimension = Hamiltonian%actual_matrix_dimension
    CALL ConstructEmptyDistributedSparseMatrix(Density, matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(WorkingHamiltonian, &
         & matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(P_k, matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(G_k, matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(G_kplusone, matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(H_k, matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(TempMat, matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(TempMat2, matrix_dimension)
    CALL ConstructEmptyDistributedSparseMatrix(Gradient, matrix_dimension)
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

    !! Initialize
    trace_value = 0.0d+0
    last_trace_value = 0.0d+0
    CALL CopyDistributedSparseMatrix(Identity, P_k)
    CALL ScaleDistributedSparseMatrix(P_k,REAL(0.5*nel,NTREAL)/matrix_dimension)

    !! Compute The Gradient
    CALL CopyDistributedSparseMatrix(Identity,TempMat)
    CALL IncrementDistributedSparseMatrix(P_k,TempMat,REAL(-1.0,NTREAL))
    CALL DistributedGemm(P_k,WorkingHamiltonian,TempMat2, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
    CALL DistributedGemm(TempMat,TempMat2,Gradient, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
    CALL ScaleDistributedSparseMatrix(Gradient,REAL(6.0,NTREAL))
    mu = Trace(Gradient)/matrix_dimension
    CALL IncrementDistributedSparseMatrix(Identity,Gradient, &
         & REAL(-1.0*mu,NTREAL))
    CALL CopyDistributedSparseMatrix(Gradient,H_k)
    CALL ScaleDistributedSparseMatrix(H_k,REAL(-1.0,NTREAL))
    CALL CopyDistributedSparseMatrix(H_k,G_k)

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0d+0
    energy_value = 0
    DO outer_counter = 1,solver_parameters%max_iterations
       IF (solver_parameters%be_verbose .AND. outer_counter .GT. 1) THEN
          CALL WriteListElement(key="Round", int_value_in=outer_counter-1)
          CALL EnterSubLog
          CALL WriteListElement(key="Convergence", float_value_in=norm_value)
          CALL WriteListElement("Energy_Value", float_value_in=energy_value)
          CALL ExitSubLog
       END IF

       !! Compute The Step Size
       b = DotDistributedSparseMatrix(H_k,Gradient)
       CALL DistributedGemm(H_k,WorkingHamiltonian,TempMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
       CALL DistributedGemm(H_k,TempMat,TempMat2, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)

       CALL CopyDistributedSparseMatrix(Identity,TempMat)
       CALL IncrementDistributedSparseMatrix(P_k,TempMat,REAL(-2.0,NTREAL))
       c = 3.0*DotDistributedSparseMatrix(TempMat,TempMat2)
       d = -2.0*DotDistributedSparseMatrix(H_k,TempMat2)

       !! Find optimal step size by solving a quadratic equation.
       root_temp = SQRT(c*c - 3.0 * b * d)
       root1 = (-1.0*root_temp - c)/(3.0*d)
       root2 = (root_temp - c)/(3.0*d)
       IF (2.0*c + 6.0*D*root1 .GT. 0) THEN
          step_size = root1
       ELSE
          step_size = root2
       END IF

       CALL IncrementDistributedSparseMatrix(H_k,P_k,REAL(step_size,NTREAL))

       !! Compute new gradient
       CALL CopyDistributedSparseMatrix(Identity,TempMat)
       CALL IncrementDistributedSparseMatrix(P_k,TempMat,REAL(-1.0,NTREAL))
       CALL DistributedGemm(P_k,WorkingHamiltonian,TempMat2, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
       CALL DistributedGemm(TempMat,TempMat2,Gradient, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
       CALL ScaleDistributedSparseMatrix(Gradient,REAL(6.0,NTREAL))
       mu = Trace(Gradient)/matrix_dimension
       CALL IncrementDistributedSparseMatrix(Identity,Gradient, &
            & REAL(-1.0*mu,NTREAL))
       CALL CopyDistributedSparseMatrix(Gradient,G_kplusone)
       CALL ScaleDistributedSparseMatrix(G_kplusone,REAL(-1.0,NTREAL))

       !! Compute conjugate direction
       gamma = DotDistributedSparseMatrix(G_kplusone,G_kplusone)
       gamma = gamma/(DotDistributedSparseMatrix(G_k,G_k))
       IF (gamma < 0.0) THEN
          gamma = 0
       END IF
       CALL ScaleDistributedSparseMatrix(H_k,gamma)
       CALL IncrementDistributedSparseMatrix(G_kplusone,H_k)
       CALL CopyDistributedSparseMatrix(G_kplusone,G_k)

       !! Energy value based convergence
       energy_value2 = energy_value
       energy_value = 2*DotDistributedSparseMatrix(P_k, WorkingHamiltonian)
       norm_value = ABS(energy_value - energy_value2)

       trace_value = Trace(P_k)

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations",int_value_in=outer_counter-1)
       CALL PrintMatrixInformation(P_k)
    END IF

    !! Undo Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(P_k, P_k, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
    END IF

    !! Compute the density matrix in the non-orthogonalized basis
    CALL DistributedGemm(InverseSquareRoot,P_k,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
    CALL DistributedGemm(TempMat,InverseSquareRoot,Density, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)

    IF (PRESENT(chemical_potential_out)) THEN
       chemical_potential_out = mu
    END IF

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    CALL DestructDistributedSparseMatrix(WorkingHamiltonian)
    CALL DestructDistributedSparseMatrix(Identity)
    CALL DestructDistributedSparseMatrix(TempMat)
    CALL DestructDistributedSparseMatrix(TempMat2)
    CALL DestructDistributedSparseMatrix(Gradient)
    CALL DestructDistributedSparseMatrix(P_k)
    CALL DestructDistributedSparseMatrix(G_k)
    CALL DestructDistributedSparseMatrix(G_kplusone)
    CALL DestructDistributedSparseMatrix(H_k)
    CALL DestructDistributedMatrixMemoryPool(pool1)
  END SUBROUTINE ConjugateGradient
END MODULE MinimizerSolversModule
