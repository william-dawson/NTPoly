!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing The Inverse of a Matrix.
MODULE InverseSolversModule
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
       & DestructSolverParameters
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
  SUBROUTINE Invert(Mat, InverseMat, solver_parameters_in)
    !> The matrix to invert.
    TYPE(Matrix_ps), INTENT(IN)  :: Mat
    !> The inverse of that matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: InverseMat
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Variables
    REAL(NTREAL) :: sigma
    TYPE(Matrix_ps) :: Temp1,Temp2,Identity
    TYPE(Matrix_ps) :: BalancedMat
    !! Temporary Variables
    INTEGER :: outer_counter
    REAL(NTREAL) :: norm_value
    TYPE(MatrixMemoryPool_p) :: pool

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Inverse Solver")
       CALL EnterSubLog
       CALL WriteHeader("Citations")
       CALL EnterSubLog
       CALL WriteListElement("palser1998canonical")
       CALL ExitSubLog
       CALL PrintParameters(solver_parameters)
    END IF

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(InverseMat, Mat)
    CALL ConstructEmptyMatrix(Temp1, Mat)
    CALL ConstructEmptyMatrix(Temp2, Mat)
    CALL ConstructEmptyMatrix(Identity, Mat)
    CALL ConstructEmptyMatrix(BalancedMat, Mat)
    CALL FillMatrixIdentity(Identity)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(Mat, BalancedMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    ELSE
       CALL CopyMatrix(Mat,BalancedMat)
    END IF

    !! Compute Sigma
    CALL MatrixSigma(BalancedMat,sigma)

    !! Create Inverse Guess
    CALL CopyMatrix(BalancedMat,InverseMat)
    CALL ScaleMatrix(InverseMat,sigma)

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0_NTREAL
    DO outer_counter = 1,solver_parameters%max_iterations
       IF (solver_parameters%be_verbose .AND. outer_counter .GT. 1) THEN
          CALL WriteListElement(key="Convergence", VALUE=norm_value)
       END IF

       CALL MatrixMultiply(InverseMat,BalancedMat,Temp1, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)

       !! Check if Converged
       CALL CopyMatrix(Identity,Temp2)
       CALL IncrementMatrix(Temp1,Temp2,-1.0_NTREAL)
       norm_value = MatrixNorm(Temp2)

       CALL DestructMatrix(Temp2)
       CALL MatrixMultiply(Temp1,InverseMat,Temp2,alpha_in=-1.0_NTREAL, &
            & threshold_in=solver_parameters%threshold,memory_pool_in=pool)

       !! Save a copy of the last inverse matrix
       CALL CopyMatrix(InverseMat,Temp1)

       CALL ScaleMatrix(InverseMat,2.0_NTREAL)

       CALL IncrementMatrix(Temp2,InverseMat, &
            & threshold_in=solver_parameters%threshold)

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations", VALUE=outer_counter-1)
       CALL PrintMatrixInformation(InverseMat)
    END IF

    !! Undo Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(InverseMat,InverseMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    !! Cleanup
    CALL DestructMatrix(Temp1)
    CALL DestructMatrix(Temp2)
    CALL DestructMatrix(BalancedMat)
    CALL DestructMatrixMemoryPool(pool)
    CALL DestructSolverParameters(solver_parameters)
  END SUBROUTINE Invert
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the inverse of a matrix using the eigendecomposition.
  SUBROUTINE DenseInvert(Mat, InverseMat, solver_parameters_in)
    !> The matrix to compute the pseudo inverse of.
    TYPE(Matrix_ps), INTENT(IN)  :: Mat
    !> The pseudoinverse of the input matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: InverseMat
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Inverse Solver")
       CALL EnterSubLog
       CALL PrintParameters(solver_parameters)
    END IF

    !! Apply
    CALL DenseMatrixFunction(Mat, InverseMat, InvertLambda, solver_parameters)

    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF
  END SUBROUTINE DenseInvert
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the pseudoinverse of a matrix.
  !> An implementation of the method of Hotelling \cite palser1998canonical.
  SUBROUTINE PseudoInverse(Mat, InverseMat, solver_parameters_in)
    !> The matrix to compute the pseudo inverse of.
    TYPE(Matrix_ps), INTENT(IN)  :: Mat
    !> The pseudoinverse of the input matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: InverseMat
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Variables
    REAL(NTREAL) :: sigma
    TYPE(Matrix_ps) :: Temp1,Temp2,Identity
    TYPE(Matrix_ps) :: BalancedMat
    !! Temporary Variables
    INTEGER :: outer_counter
    REAL(NTREAL) :: norm_value
    TYPE(MatrixMemoryPool_p) :: pool

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Inverse Solver")
       CALL EnterSubLog
       CALL WriteHeader("Citations")
       CALL EnterSubLog
       CALL WriteListElement("palser1998canonical")
       CALL ExitSubLog
       CALL PrintParameters(solver_parameters)
    END IF

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(InverseMat, Mat)
    CALL ConstructEmptyMatrix(Temp1, Mat)
    CALL ConstructEmptyMatrix(Temp2, Mat)
    CALL ConstructEmptyMatrix(Identity, Mat)
    CALL ConstructEmptyMatrix(BalancedMat, Mat)
    CALL FillMatrixIdentity(Identity)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(Mat, BalancedMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    ELSE
       CALL CopyMatrix(Mat,BalancedMat)
    END IF

    !! Compute Sigma
    CALL MatrixSigma(BalancedMat,sigma)

    !! Create Inverse Guess
    CALL CopyMatrix(BalancedMat,InverseMat)
    CALL ScaleMatrix(InverseMat,sigma)

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0_NTREAL
    DO outer_counter = 1,solver_parameters%max_iterations
       IF (solver_parameters%be_verbose .AND. outer_counter .GT. 1) THEN
          CALL WriteListElement(key="Convergence", VALUE=norm_value)
       END IF

       CALL MatrixMultiply(InverseMat,BalancedMat,Temp1, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       CALL MatrixMultiply(Temp1,InverseMat,Temp2,alpha_in=-1.0_NTREAL, &
            & threshold_in=solver_parameters%threshold,memory_pool_in=pool)

       !! Save a copy of the last inverse matrix
       CALL CopyMatrix(InverseMat,Temp1)

       CALL ScaleMatrix(InverseMat,2.0_NTREAL)
       CALL IncrementMatrix(Temp2,InverseMat, &
            & threshold_in=solver_parameters%threshold)

       !! Check if Converged
       CALL IncrementMatrix(InverseMat,Temp1,-1.0_NTREAL)
       norm_value = MatrixNorm(Temp1)

       !! Sometimes the first few values don't change so much, so that's why
       !! I added the outer counter check
       IF (norm_value .LE. solver_parameters%converge_diff .AND. &
            & outer_counter .GT. 3) THEN
          EXIT
       END IF
    END DO
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations", VALUE=outer_counter-1)
       CALL PrintMatrixInformation(InverseMat)
    END IF

    !! Undo Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(InverseMat,InverseMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    !! Cleanup
    CALL DestructMatrix(Temp1)
    CALL DestructMatrix(Temp2)
    CALL DestructMatrix(BalancedMat)
    CALL DestructMatrixMemoryPool(pool)
    CALL DestructSolverParameters(solver_parameters)
  END SUBROUTINE PseudoInverse
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Prototypical inversion for mapping. 
  SUBROUTINE InvertLambda(index, val)
    !> The index of the eigenvalue
    INTEGER, INTENT(IN) :: index
    !> The actual value of an element.
    REAL(KIND=NTREAL), INTENT(INOUT) :: val

    val = 1.0 / val
  END SUBROUTINE InvertLambda
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE InverseSolversModule
