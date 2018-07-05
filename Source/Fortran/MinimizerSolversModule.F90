!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Solving Systems Quantum Chemistry Systems Using Minimization.
MODULE MinimizerSolversModule
  USE DataTypesModule
  USE LoadBalancerModule
  USE LoggingModule
  USE PMatrixMemoryPoolModule
  USE PSMatrixAlgebraModule
  USE PSMatrixModule
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters
  USE TimerModule
  USE MPI
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
  !! @param[out] chemical_potential_out the chemical potential (optional).
  !! @param[in] solver_parameters_in parameters for the solver.
  SUBROUTINE ConjugateGradient(Hamiltonian, InverseSquareRoot, nel, Density, &
       & chemical_potential_out, solver_parameters_in)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN)  :: Hamiltonian
    TYPE(Matrix_ps), INTENT(IN)  :: InverseSquareRoot
    INTEGER, INTENT(IN) :: nel
    TYPE(Matrix_ps), INTENT(INOUT) :: Density
    REAL(NTREAL), INTENT(out), OPTIONAL :: chemical_potential_out
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: &
         & solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Variables
    REAL(NTREAL) :: trace_value
    REAL(NTREAL) :: last_trace_value
    REAL(NTREAL) :: norm_value
    TYPE(Matrix_ps) :: WorkingHamiltonian
    TYPE(Matrix_ps) :: P_k
    TYPE(Matrix_ps) :: Gradient
    TYPE(Matrix_ps) :: G_k, G_kplusone
    TYPE(Matrix_ps) :: H_k
    TYPE(Matrix_ps) :: TempMat, TempMat2
    TYPE(Matrix_ps) :: Identity
    REAL(NTREAL) :: mu
    REAL(NTREAL) :: gamma, gamma_d
    REAL(NTREAL) :: step_size
    REAL(NTREAL) :: b, c, d
    REAL(NTREAL) :: root1, root2, root_temp
    REAL(NTREAL) :: energy_value, energy_value2
    !! Temporary Variables
    TYPE(MatrixMemoryPool_p) :: pool1
    INTEGER :: outer_counter
    INTEGER :: matrix_dimension

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Density Matrix Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="CG")
       CALL WriteCitation("millam1997linear challacombe1999simplified")
       CALL PrintParameters(solver_parameters)
    END IF

    !! Construct All The Necessary Matrices
    matrix_dimension = Hamiltonian%actual_matrix_dimension
    CALL ConstructEmptyMatrix(Density, Hamiltonian)
    CALL ConstructEmptyMatrix(WorkingHamiltonian, Hamiltonian)
    CALL ConstructEmptyMatrix(P_k, Hamiltonian)
    CALL ConstructEmptyMatrix(G_k, Hamiltonian)
    CALL ConstructEmptyMatrix(G_kplusone, Hamiltonian)
    CALL ConstructEmptyMatrix(H_k, Hamiltonian)
    CALL ConstructEmptyMatrix(TempMat, Hamiltonian)
    CALL ConstructEmptyMatrix(TempMat2, Hamiltonian)
    CALL ConstructEmptyMatrix(Gradient, Hamiltonian)
    CALL ConstructEmptyMatrix(Identity, Hamiltonian)
    CALL FillMatrixIdentity(Identity)

    !! Compute the working hamiltonian.
    CALL MatrixMultiply(InverseSquareRoot,Hamiltonian,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
    CALL MatrixMultiply(TempMat,InverseSquareRoot,WorkingHamiltonian, &
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
    CALL CopyMatrix(Identity, P_k)
    CALL ScaleMatrix(P_k,REAL(0.5*nel,NTREAL)/matrix_dimension)

    !! Compute The Gradient
    CALL CopyMatrix(Identity,TempMat)
    CALL IncrementMatrix(P_k,TempMat,REAL(-1.0,NTREAL))
    CALL MatrixMultiply(P_k,WorkingHamiltonian,TempMat2, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
    CALL MatrixMultiply(TempMat,TempMat2,Gradient, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
    CALL ScaleMatrix(Gradient,REAL(6.0,NTREAL))
    CALL MatrixTrace(Gradient, mu)
    mu = mu/matrix_dimension
    CALL IncrementMatrix(Identity,Gradient, &
         & REAL(-1.0*mu,NTREAL))
    CALL CopyMatrix(Gradient,H_k)
    CALL ScaleMatrix(H_k,REAL(-1.0,NTREAL))
    CALL CopyMatrix(H_k,G_k)

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
       CALL DotMatrix(H_k, Gradient, b)
       CALL MatrixMultiply(H_k,WorkingHamiltonian,TempMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
       CALL MatrixMultiply(H_k,TempMat,TempMat2, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)

       CALL CopyMatrix(Identity,TempMat)
       CALL IncrementMatrix(P_k,TempMat,REAL(-2.0,NTREAL))
       CALL DotMatrix(TempMat, TempMat2, c)
       c = 3.0*c
       CALL DotMatrix(H_k, TempMat2, d)
       d = -2.0*d

       !! Find optimal step size by solving a quadratic equation.
       root_temp = SQRT(c*c - 3.0 * b * d)
       root1 = (-1.0*root_temp - c)/(3.0*d)
       root2 = (root_temp - c)/(3.0*d)
       IF (2.0*c + 6.0*D*root1 .GT. 0) THEN
          step_size = root1
       ELSE
          step_size = root2
       END IF

       CALL IncrementMatrix(H_k,P_k,REAL(step_size,NTREAL))

       !! Compute new gradient
       CALL CopyMatrix(Identity,TempMat)
       CALL IncrementMatrix(P_k,TempMat,REAL(-1.0,NTREAL))
       CALL MatrixMultiply(P_k,WorkingHamiltonian,TempMat2, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
       CALL MatrixMultiply(TempMat,TempMat2,Gradient, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
       CALL ScaleMatrix(Gradient,REAL(6.0,NTREAL))
       CALL MatrixTrace(Gradient, mu)
       mu = mu/matrix_dimension
       CALL IncrementMatrix(Identity,Gradient, &
            & REAL(-1.0*mu,NTREAL))
       CALL CopyMatrix(Gradient,G_kplusone)
       CALL ScaleMatrix(G_kplusone,REAL(-1.0,NTREAL))

       !! Compute conjugate direction
       CALL DotMatrix(G_kplusone,G_kplusone,gamma)
       CALL DotMatrix(G_k, G_k, gamma_d)
       gamma = gamma/gamma_d
       IF (gamma < 0.0) THEN
          gamma = 0
       END IF
       CALL ScaleMatrix(H_k,gamma)
       CALL IncrementMatrix(G_kplusone,H_k)
       CALL CopyMatrix(G_kplusone,G_k)

       !! Energy value based convergence
       energy_value2 = energy_value
       CALL DotMatrix(P_k, WorkingHamiltonian, energy_value)
       energy_value = 2*energy_value
       norm_value = ABS(energy_value - energy_value2)

       CALL MatrixTrace(P_k, trace_value)

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
    CALL MatrixMultiply(InverseSquareRoot,P_k,TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
    CALL MatrixMultiply(TempMat,InverseSquareRoot,Density, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)

    IF (PRESENT(chemical_potential_out)) THEN
       chemical_potential_out = mu
    END IF

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    CALL DestructMatrix(WorkingHamiltonian)
    CALL DestructMatrix(Identity)
    CALL DestructMatrix(TempMat)
    CALL DestructMatrix(TempMat2)
    CALL DestructMatrix(Gradient)
    CALL DestructMatrix(P_k)
    CALL DestructMatrix(G_k)
    CALL DestructMatrix(G_kplusone)
    CALL DestructMatrix(H_k)
    CALL DestructMatrixMemoryPool(pool1)
  END SUBROUTINE ConjugateGradient
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE MinimizerSolversModule
