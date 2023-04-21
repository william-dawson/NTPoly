!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing General Matrix Roots.
MODULE RootSolversModule
  USE DataTypesModule, ONLY : NTREAL
  USE EigenBoundsModule, ONLY : GershgorinBounds
  USE InverseSolversModule, ONLY : Invert
  USE LoadBalancerModule, ONLY : PermuteMatrix, UndoPermuteMatrix
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteListElement, &
       & WriteHeader, WriteElement
  USE PolynomialSolversModule, ONLY : Polynomial_t, ConstructPolynomial, &
       & DestructPolynomial, FactorizedCompute, SetCoefficient
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply, MatrixNorm, &
       & IncrementMatrix, ScaleMatrix
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, CopyMatrix, &
       & DestructMatrix, FillMatrixIdentity, PrintMatrixInformation
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters, &
       & DestructSolverParameters, ConstructSolverParameters, &
       & CopySolverParameters
  USE SquareRootSolversModule, ONLY : SquareRoot, InverseSquareRoot
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Solvers
  PUBLIC :: ComputeRoot
  PUBLIC :: ComputeInverseRoot
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a general matrix root.
  RECURSIVE SUBROUTINE ComputeRoot(InputMat, OutputMat, root, &
       & solver_parameters_in)
    !> The input matrix
    TYPE(Matrix_ps), INTENT(IN)  :: InputMat
    !> OutputMat = InputMat^1/root.
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> Which root to compute.
    INTEGER, INTENT(IN) :: root
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(SolverParameters_t) :: params
    !! Local Variables
    TYPE(Matrix_ps) :: TempMat

    !! Handle The Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF

    IF (params%be_verbose) THEN
       CALL WriteHeader("Root Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Root", VALUE=root)
       CALL PrintParameters(params)
    END IF

    !! Handle base cases, or call to general implementation.
    IF (root .EQ. 1) THEN
       CALL CopyMatrix(InputMat, OutputMat)
    ELSE IF (root .EQ. 2) THEN
       CALL SquareRoot(InputMat, OutputMat, params)
    ELSE IF (root .EQ. 3) THEN
       CALL MatrixMultiply(InputMat, InputMat, TempMat, &
            & threshold_in=params%threshold)
       CALL ComputeRootImplementation(TempMat, OutputMat, 6, &
            & params)
    ELSE IF (root .EQ. 4) THEN
       CALL SquareRoot(InputMat, TempMat, params)
       CALL SquareRoot(TempMat, OutputMat, params)
       CALL DestructMatrix(TempMat)
    ELSE
       CALL ComputeRootImplementation(InputMat, OutputMat, root, &
            & params)
    END IF

    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF

    CALL DestructSolverParameters(params)

  END SUBROUTINE ComputeRoot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Actual implementation of computing a general matrix root.
  SUBROUTINE ComputeRootImplementation(InputMat, OutputMat, root, params)
    !> The input matrix
    TYPE(Matrix_ps), INTENT(IN)  :: InputMat
    !> OutputMat = InputMat^1/root.
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> Which root to compute.
    INTEGER, INTENT(IN) :: root
    !> Parameters for the solver
    TYPE(SolverParameters_t), INTENT(IN) :: params
    !! Local Variables
    TYPE(Matrix_ps) :: RaisedMat
    TYPE(Matrix_ps) :: TempMat
    TYPE(Polynomial_t) :: power_poly
    INTEGER :: II

    !! We will use the formula A^(1/x) = A*A^(1/x - 1)
    !! So first, we raise to the root-1 power
    CALL ConstructPolynomial(power_poly, root)
    DO II = 1, root-1
       CALL SetCoefficient(power_poly, II, REAL(0.0,NTREAL))
    END DO
    CALL SetCoefficient(power_poly, root, REAL(1.0,NTREAL))
    CALL FactorizedCompute(InputMat, RaisedMat, power_poly, params)
    CALL DestructPolynomial(power_poly)

    !! Now compute the inverse pth root
    CALL ComputeInverseRoot(RaisedMat, TempMat, root, params)

    !! Multiply by the original matrix
    CALL MatrixMultiply(InputMat, TempMat, OutputMat, &
         & threshold_in=params%threshold)

    !! Cleanup
    CALL DestructMatrix(RaisedMat)
    CALL DestructMatrix(TempMat)
  END SUBROUTINE ComputeRootImplementation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a general inverse matrix root.
  RECURSIVE SUBROUTINE ComputeInverseRoot(InputMat, OutputMat, root, &
       & solver_parameters_in)
    !> The input matrix
    TYPE(Matrix_ps), INTENT(IN)  :: InputMat
    !> OutputMat = InputMat^-1/root.
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> Which root to compute.
    INTEGER, INTENT(IN) :: root
    !> Parameters for the solver.
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Solver Parameters
    TYPE(SolverParameters_t) :: params
    !! Local Variables
    TYPE(Matrix_ps) :: TempMat

    !! Handle The Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       CALL CopySolverParameters(solver_parameters_in, params)
    ELSE
       CALL ConstructSolverParameters(params)
    END IF

    IF (params%be_verbose) THEN
       CALL WriteHeader("Inverse Root Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Root", VALUE=root)
       CALL PrintParameters(params)
    END IF

    !! Handle base cases, or call to general implementation.
    IF (root .EQ. 1) THEN
       CALL Invert(InputMat, OutputMat, params)
    ELSE IF (root .EQ. 2) THEN
       CALL InverseSquareRoot(InputMat, OutputMat, params)
    ELSE IF (root .EQ. 3) THEN
       CALL ComputeRoot(InputMat,TempMat,3, params)
       CALL Invert(TempMat, OutputMat, params)
    ELSE IF (root .EQ. 4) THEN
       CALL SquareRoot(InputMat, TempMat, params)
       CALL InverseSquareRoot(TempMat, OutputMat, params)
       CALL DestructMatrix(TempMat)
    ELSE
       CALL ComputeInverseRootImplemention(InputMat, OutputMat, root, params)
    END IF

    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF

    CALL DestructSolverParameters(params)
  END SUBROUTINE ComputeInverseRoot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a general inverse matrix root for root > 4.
  SUBROUTINE ComputeInverseRootImplemention(InputMat, OutputMat, root, params)
    !> Matrix to compute the root of.
    TYPE(Matrix_ps), INTENT(IN)  :: InputMat
    !> The inverse nth root of that matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: OutputMat
    !> Which inverse root to compute.
    INTEGER, INTENT(IN) :: root
    !> Parameters for the solver.
    TYPE(SolverParameters_t), INTENT(IN) :: params
    !! Constants.
    REAL(NTREAL), PARAMETER :: NEGATIVE_ONE = -1.0
    !! Local Matrices
    TYPE(Matrix_ps) :: SqrtMat, FthrtMat
    TYPE(Matrix_ps) :: IdentityMat
    TYPE(Matrix_ps) :: Mk
    TYPE(Matrix_ps) :: IntermediateMat
    TYPE(Matrix_ps) :: IntermediateMatP
    TYPE(Matrix_ps) :: Temp
    !! Local Variables
    INTEGER :: target_root
    REAL(NTREAL) :: e_min
    REAL(NTREAL) :: e_max
    REAL(NTREAL) :: scaling_factor
    REAL(NTREAL) :: norm_value
    !! Temporary Variables
    INTEGER :: II
    INTEGER :: JJ
    TYPE(MatrixMemoryPool_p) :: pool

    IF (params%be_verbose) THEN
       CALL WriteHeader("Root Solver")
       CALL EnterSubLog
       CALL WriteHeader("Citations")
       CALL EnterSubLog
       CALL WriteListElement("nicholas2008functions")
       CALL ExitSubLog
       CALL PrintParameters(params)
    END IF

    !! Compute The Scaling Factor
    CALL GershgorinBounds(InputMat, e_min, e_max)
    scaling_factor = e_max/SQRT(2.0)**(1.0/root)

    !! Compute the target root (adjust for the fact that we just took the
    !! fourth root.
    target_root = 0
    IF (MOD(root,4) .EQ. 0) THEN
       target_root = root/4
    ELSE IF (MOD(root,4) .EQ. 1 .OR. MOD(root,4) .EQ. 3) THEN
       target_root = root
    ELSE
       target_root = (root-2)/2 + 1
    END IF

    !! Initialize
    !! Fourth Root Matrix
    CALL SquareRoot(InputMat, SqrtMat, params)
    CALL SquareRoot(SqrtMat, FthrtMat, params)
    CALL DestructMatrix(SqrtMat)

    !! Setup the Matrices
    CALL ConstructEmptyMatrix(IdentityMat, InputMat)
    CALL FillMatrixIdentity(IdentityMat)

    !! Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL PermuteMatrix(FthrtMat, FthrtMat, &
            & params%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(IdentityMat, IdentityMat, &
            & params%BalancePermutation, memorypool_in=pool)
    END IF

    CALL CopyMatrix(IdentityMat, OutputMat)
    CALL ScaleMatrix(OutputMat, 1.0/scaling_factor)

    CALL CopyMatrix(FthrtMat, Mk)
    CALL ScaleMatrix(Mk, 1.0/(scaling_factor**target_root))
    CALL DestructMatrix(FthrtMat)

    CALL ConstructEmptyMatrix(IntermediateMat, InputMat)
    CALL ConstructEmptyMatrix(IntermediateMatP, InputMat)
    CALL ConstructEmptyMatrix(Temp, InputMat)

    IF (params%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    II = 1
    norm_value = params%converge_diff + 1.0_NTREAL
    DO II = 1, params%max_iterations
       IF (params%be_verbose .AND. II .GT. 1) THEN
          CALL WriteListElement(key="Convergence", VALUE=norm_value)
       END IF

       CALL CopyMatrix(IdentityMat, IntermediateMat)
       CALL ScaleMatrix(IntermediateMat, &
            & REAL(target_root+1,NTREAL))
       CALL IncrementMatrix(Mk, IntermediateMat, &
            & alpha_in=NEGATIVE_ONE)
       CALL ScaleMatrix(IntermediateMat, &
            & REAL(1.0,NTREAL)/target_root)

       CALL MatrixMultiply(OutputMat, IntermediateMat, Temp, &
            & threshold_in=params%threshold, memory_pool_in=pool)
       CALL CopyMatrix(Temp, OutputMat)

       CALL CopyMatrix(IntermediateMat, IntermediateMatP)
       DO JJ = 1, target_root-1
          CALL MatrixMultiply(IntermediateMat, IntermediateMatP, Temp, &
               & threshold_in=params%threshold, memory_pool_in=pool)
          CALL CopyMatrix(Temp, IntermediateMatP)
       END DO

       CALL MatrixMultiply(IntermediateMatP, Mk, Temp, &
            & threshold_in=params%threshold, memory_pool_in=pool)
       CALL CopyMatrix(Temp, Mk)

       CALL IncrementMatrix(IdentityMat, Temp, &
            & alpha_in=NEGATIVE_ONE)
       norm_value = MatrixNorm(Temp)

       IF (norm_value .LE. params%converge_diff) THEN
          EXIT
       END IF
    END DO
    IF (params%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations", VALUE=II-1)
       CALL PrintMatrixInformation(OutputMat)
    END IF

    IF (MOD(root,4) .EQ. 1 .OR. MOD(root,4) .EQ. 3) THEN
       CALL MatrixMultiply(OutputMat, OutputMat, Temp, &
            & threshold_in=params%threshold, memory_pool_in=pool)
       CALL MatrixMultiply(Temp, Temp, OutputMat, &
            & threshold_in=params%threshold, memory_pool_in=pool)
    ELSE IF (MOD(root,4) .NE. 0) THEN
       CALL MatrixMultiply(OutputMat, OutputMat, Temp, &
            & threshold_in=params%threshold, memory_pool_in=pool)
       CALL CopyMatrix(Temp, OutputMat)
    END IF

    !! Undo Load Balancing Step
    IF (params%do_load_balancing) THEN
       CALL UndoPermuteMatrix(OutputMat, OutputMat, &
            & params%BalancePermutation, memorypool_in=pool)
    END IF

    !! Cleanup
    IF (params%be_verbose) THEN
       CALL ExitSubLog
    END IF
    CALL DestructMatrix(IdentityMat)
    CALL DestructMatrix(Mk)
    CALL DestructMatrix(IntermediateMat)
    CALL DestructMatrix(IntermediateMatP)
    CALL DestructMatrix(Temp)
    CALL DestructMatrixMemoryPool(pool)
  END SUBROUTINE ComputeInverseRootImplemention
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE RootSolversModule
