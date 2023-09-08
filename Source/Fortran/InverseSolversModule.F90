!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing The Inverse of a Matrix.
MODULE InverseSolversModule
  USE ConvergenceMonitor, ONLY : ConstructMonitor, CheckConverged, AppendValue
  USE DataTypesModule, ONLY : NTREAL
  USE EigenSolversModule, ONLY : DenseMatrixFunction
  USE LoadBalancerModule, ONLY : PermuteMatrix, UndoPermuteMatrix
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteHeader, &
       & WriteElement, WriteListElement
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply, IncrementMatrix, &
       & MatrixNorm, MatrixSigma, ScaleMatrix
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, CopyMatrix, &
       & DestructMatrix, FillMatrixIdentity, PrintMatrixInformation
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters, &
       & DestructSolverParameters, CopySolverParameters, &
       & ConstructSolverParameters
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Solvers
  PUBLIC :: Invert
  PUBLIC :: DenseInvert
  PUBLIC :: PseudoInverse
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the inverse of a matrix.
  !> An implementation of the method of Hotelling \cite palser1998canonical.
  SUBROUTINE Invert(InputMat, OutputMat, solver_parameters_in)
    !> The matrix to invert.
    TYPE(Matrix_ps), INTENT(IN)  :: InputMat
    !> The inverse of that matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: params
    !! Local Variables
    REAL(NTREAL) :: sigma
    TYPE(Matrix_ps) :: Temp1,Temp2,Identity
    TYPE(Matrix_ps) :: BalancedMat
    !! Temporary Variables
    INTEGER :: II
    REAL(NTREAL) :: norm_value
    TYPE(MatrixMemoryPool_p) :: pool

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF
    CALL ConstructMonitor(params%monitor, &
         & automatic_in = params%monitor_convergence, &
         & tight_cutoff_in=params%converge_diff)

    IF (params%be_verbose) THEN
       CALL WriteHeader("Inverse Solver")
       CALL EnterSubLog
       CALL WriteHeader("Citations")
       CALL EnterSubLog
       CALL WriteListElement("palser1998canonical")
       CALL ExitSubLog
       CALL PrintParameters(params)
    END IF

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(OutputMat, InputMat)
    CALL ConstructEmptyMatrix(Temp1, InputMat)
    CALL ConstructEmptyMatrix(Temp2, InputMat)
    CALL ConstructEmptyMatrix(Identity, InputMat)
    CALL ConstructEmptyMatrix(BalancedMat, InputMat)
    CALL FillMatrixIdentity(Identity)

    !! Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL PermuteMatrix(Identity, Identity, &
            & params%BalancePermutation, memorypool_in = pool)
       CALL PermuteMatrix(InputMat, BalancedMat, &
            & params%BalancePermutation, memorypool_in = pool)
    ELSE
       CALL CopyMatrix(InputMat, BalancedMat)
    END IF

    !! Compute Sigma
    CALL MatrixSigma(BalancedMat, sigma)

    !! Create Inverse Guess
    CALL CopyMatrix(BalancedMat, OutputMat)
    CALL ScaleMatrix(OutputMat, sigma)

    !! Iterate
    IF (params%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    II = 1
    norm_value = params%converge_diff + 1.0_NTREAL
    DO II = 1, params%max_iterations
       IF (params%be_verbose .AND. II .GT. 1) THEN
          CALL WriteListElement(key = "Convergence", VALUE = norm_value)
       END IF

       CALL MatrixMultiply(OutputMat, BalancedMat, Temp1, &
            & threshold_in = params%threshold, memory_pool_in = pool)

       !! Check if Converged
       CALL CopyMatrix(Identity, Temp2)
       CALL IncrementMatrix(Temp1, Temp2, -1.0_NTREAL)
       norm_value = MatrixNorm(Temp2)

       CALL DestructMatrix(Temp2)
       CALL MatrixMultiply(Temp1, OutputMat, Temp2, alpha_in = -1.0_NTREAL, &
            & threshold_in = params%threshold, memory_pool_in = pool)

       !! Save a copy of the last inverse matrix
       CALL CopyMatrix(OutputMat, Temp1)

       CALL ScaleMatrix(OutputMat, 2.0_NTREAL)

       CALL IncrementMatrix(Temp2, OutputMat, &
            & threshold_in = params%threshold)

       CALL AppendValue(params%monitor, norm_value)
       IF (CheckConverged(params%monitor, params%be_verbose)) EXIT
    END DO
    IF (params%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key = "Total Iterations", VALUE = II-1)
       CALL PrintMatrixInformation(OutputMat)
    END IF

    !! Undo Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL UndoPermuteMatrix(OutputMat, OutputMat, &
            & params%BalancePermutation, memorypool_in = pool)
    END IF

    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF

    !! Cleanup
    CALL DestructMatrix(Temp1)
    CALL DestructMatrix(Temp2)
    CALL DestructMatrix(BalancedMat)
    CALL DestructMatrixMemoryPool(pool)
    CALL DestructSolverParameters(params)
  END SUBROUTINE Invert
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the inverse of a matrix using the eigendecomposition.
  SUBROUTINE DenseInvert(InputMat, OutputMat, solver_parameters_in)
    !> The matrix to compute the pseudo inverse of.
    TYPE(Matrix_ps), INTENT(IN)  :: InputMat
    !> The pseudoinverse of the input matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: params

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF

    IF (params%be_verbose) THEN
       CALL WriteHeader("Inverse Solver")
       CALL EnterSubLog
    END IF

    !! Apply
    CALL DenseMatrixFunction(InputMat, OutputMat, InvertLambda, params)

    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF

    !! Cleanup
    CALL DestructSolverParameters(params)
  END SUBROUTINE DenseInvert
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the pseudoinverse of a matrix.
  !> An implementation of the method of Hotelling \cite palser1998canonical.
  SUBROUTINE PseudoInverse(InputMat, OutputMat, solver_parameters_in)
    !> The matrix to compute the pseudo inverse of.
    TYPE(Matrix_ps), INTENT(IN)  :: InputMat
    !> The pseudoinverse of the input matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: params
    !! Local Variables
    REAL(NTREAL) :: sigma
    TYPE(Matrix_ps) :: Temp1,Temp2,Identity
    TYPE(Matrix_ps) :: BalancedMat
    !! Temporary Variables
    INTEGER :: II
    REAL(NTREAL) :: norm_value
    TYPE(MatrixMemoryPool_p) :: pool

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF
    CALL ConstructMonitor(params%monitor, &
         & automatic_in = params%monitor_convergence, &
         & tight_cutoff_in=params%converge_diff)

    IF (params%be_verbose) THEN
       CALL WriteHeader("Inverse Solver")
       CALL EnterSubLog
       CALL WriteHeader("Citations")
       CALL EnterSubLog
       CALL WriteListElement("palser1998canonical")
       CALL ExitSubLog
       CALL PrintParameters(params)
    END IF

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(OutputMat, InputMat)
    CALL ConstructEmptyMatrix(Temp1, InputMat)
    CALL ConstructEmptyMatrix(Temp2, InputMat)
    CALL ConstructEmptyMatrix(Identity, InputMat)
    CALL ConstructEmptyMatrix(BalancedMat, InputMat)
    CALL FillMatrixIdentity(Identity)

    !! Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL PermuteMatrix(Identity, Identity, &
            & params%BalancePermutation, memorypool_in = pool)
       CALL PermuteMatrix(InputMat, BalancedMat, &
            & params%BalancePermutation, memorypool_in = pool)
    ELSE
       CALL CopyMatrix(InputMat, BalancedMat)
    END IF

    !! Compute Sigma
    CALL MatrixSigma(BalancedMat, sigma)

    !! Create Inverse Guess
    CALL CopyMatrix(BalancedMat, OutputMat)
    CALL ScaleMatrix(OutputMat, sigma)

    !! Iterate
    IF (params%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    II = 1
    norm_value = params%converge_diff + 1.0_NTREAL
    DO II = 1,params%max_iterations
       CALL MatrixMultiply(OutputMat, BalancedMat, Temp1, &
            & threshold_in = params%threshold, memory_pool_in = pool)
       CALL MatrixMultiply(Temp1, OutputMat, Temp2, alpha_in = -1.0_NTREAL, &
            & threshold_in = params%threshold, memory_pool_in = pool)

       !! Save a copy of the last inverse matrix
       CALL CopyMatrix(OutputMat, Temp1)

       CALL ScaleMatrix(OutputMat, 2.0_NTREAL)
       CALL IncrementMatrix(Temp2, OutputMat, &
            & threshold_in = params%threshold)

       !! Check if Converged
       CALL IncrementMatrix(OutputMat, Temp1, -1.0_NTREAL)
       norm_value = MatrixNorm(Temp1)

       CALL AppendValue(params%monitor, norm_value)
       IF (CheckConverged(params%monitor, params%be_verbose)) EXIT
    END DO
    IF (params%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key = "Total Iterations", VALUE = II - 1)
       CALL PrintMatrixInformation(OutputMat)
    END IF

    !! Undo Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL UndoPermuteMatrix(OutputMat, OutputMat, &
            & params%BalancePermutation, memorypool_in = pool)
    END IF

    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF

    !! Cleanup
    CALL DestructMatrix(Temp1)
    CALL DestructMatrix(Temp2)
    CALL DestructMatrix(BalancedMat)
    CALL DestructMatrixMemoryPool(pool)
    CALL DestructSolverParameters(params)
  END SUBROUTINE PseudoInverse
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Prototypical inversion for mapping. 
  FUNCTION InvertLambda(val) RESULT(outval)
    REAL(KIND = NTREAL), INTENT(IN) :: val
    REAL(KIND = NTREAL) :: outval

    outval = 1.0 / val
  END FUNCTION InvertLambda
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE InverseSolversModule
