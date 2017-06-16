!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing General Matrix Roots.
MODULE RootSolversModule
  USE DataTypesModule
#ifdef NOBLOCK
  USE DistributedSparseMatrixModule
#else
  USE DistributedBlockedSparseMatrixModule
#endif
  USE DistributedMatrixMemoryPoolModule
  USE FixedSolversModule
  USE InverseSolversModule
  USE IterativeSolversModule
  USE LoadBalancerModule
  USE LoggingModule
  USE ProcessGridModule
  USE PolynomialSolversModule
  USE SquareRootSolversModule
  USE TimerModule
  USE mpi
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Solvers
  PUBLIC :: ComputeRoot
  PUBLIC :: ComputeInverseRoot
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a general matrix root.
  !! @param[in] InputMat the input matrix
  !! @param[out] OutputMat = InputMat^1/root.
  !! @param[in] root which root to compute.
  !! @param[in] solver_parameters_in parameters for the solver (optional).
  SUBROUTINE ComputeRoot(InputMat, OutputMat, root, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in)  :: InputMat
    TYPE(DistributedSparseMatrix), INTENT(inout) :: OutputMat
    INTEGER, INTENT(in) :: root
    TYPE(IterativeSolverParameters), INTENT(in), OPTIONAL :: &
         & solver_parameters_in
    !! Handling Solver Parameters
    TYPE(IterativeSolverParameters) :: solver_parameters
    !! Local Variables
    TYPE(DistributedSparseMatrix) :: TempMat

    !! Handle The Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters()
    END IF

    IF (solver_parameters%be_verbose .AND. IsRoot()) THEN
       CALL WriteHeader("Root Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Root", int_value_in=root)
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    !! Handle base cases, or call to general implementation.
    IF (root .EQ. 1) THEN
       CALL CopyDistributedSparseMatrix(InputMat, OutputMat)
    ELSE IF (root .EQ. 2) THEN
       CALL SquareRoot(InputMat, OutputMat, solver_parameters)
    ELSE IF (root .EQ. 3) THEN
       WRITE(*,*) "NOT IMPLEMENTED YET!"
    ELSE IF (root .EQ. 4) THEN
       CALL SquareRoot(InputMat, TempMat, solver_parameters)
       CALL SquareRoot(TempMat, OutputMat, solver_parameters)
       CALL DestructDistributedSparseMatrix(TempMat)
    ELSE
       CALL ComputeRootImplementation(InputMat, OutputMat, root, &
            & solver_parameters)
    END IF

    IF (solver_parameters%be_verbose .AND. IsRoot()) THEN
       CALL ExitSubLog
    END IF

  END SUBROUTINE ComputeRoot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a general matrix root.
  !! @param[in] InputMat the input matrix
  !! @param[out] OutputMat = InputMat^1/root.
  !! @param[in] root which root to compute.
  !! @param[in] solver_parameters_in parameters for the solver (optional).
  SUBROUTINE ComputeRootImplementation(InputMat, OutputMat, root, &
       & solver_parameters)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in)  :: InputMat
    TYPE(DistributedSparseMatrix), INTENT(inout) :: OutputMat
    INTEGER, INTENT(in) :: root
    TYPE(IterativeSolverParameters), INTENT(in) :: solver_parameters
    !! Handling Solver Parameters
    TYPE(FixedSolverParameters) :: fixed_parameters
    !! Local Variables
    TYPE(DistributedSparseMatrix) :: RaisedMat
    TYPE(DistributedSparseMatrix) :: TempMat
    TYPE(Polynomial_t) :: power_poly
    INTEGER :: counter

    !! Set up the solver parameters
    fixed_parameters%threshold = solver_parameters%threshold
    fixed_parameters%be_verbose = solver_parameters%be_verbose
    fixed_parameters%do_load_balancing = solver_parameters%do_load_balancing
    fixed_parameters%BalancePermutation = solver_parameters%BalancePermutation

    !! We will use the formula A^(1/x) = A*A^(1/x - 1)
    !! So first, we raise to the root-1 power
    CALL ConstructPolynomial(power_poly,root)
    DO counter=1,root-1
       CALL SetCoefficient(power_poly,counter,REAL(0.0,NTREAL))
    END DO
    CALL SetCoefficient(power_poly,root,REAL(1.0,NTREAL))
    CALL PatersonStockmeyerCompute(InputMat, RaisedMat, power_poly, &
         & fixed_parameters)
    CALL DestructPolynomial(power_poly)

    !! Now compute the inverse pth root
    CALL ComputeInverseRoot(RaisedMat, TempMat, root, solver_parameters)

    !! Multiply by the original matrix
    CALL DistributedGemm(InputMat, TempMat, OutputMat, &
         & threshold_in=solver_parameters%threshold)

    !! Cleanup
    CALL DestructDistributedSparseMatrix(RaisedMat)
    CALL DestructDistributedSparseMatrix(TempMat)
  END SUBROUTINE ComputeRootImplementation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a general inverse matrix root.
  !! @param[in] InputMat the input matrix
  !! @param[out] OutputMat = InputMat^-1/root.
  !! @param[in] root which root to compute.
  !! @param[in] solver_parameters_in parameters for the solver (optional).
  !! @todo the inverse cube root.
  SUBROUTINE ComputeInverseRoot(InputMat, OutputMat, root, solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in)  :: InputMat
    TYPE(DistributedSparseMatrix), INTENT(inout) :: OutputMat
    INTEGER, INTENT(in) :: root
    TYPE(IterativeSolverParameters), INTENT(in), OPTIONAL :: &
         & solver_parameters_in
    !! Handling Solver Parameters
    TYPE(IterativeSolverParameters) :: solver_parameters
    !! Local Variables
    TYPE(DistributedSparseMatrix) :: TempMat

    !! Handle The Optional Parameters
    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters()
    END IF

    IF (solver_parameters%be_verbose .AND. IsRoot()) THEN
       CALL WriteHeader("Inverse Root Solver")
       CALL EnterSubLog
       CALL WriteElement(key="Root", int_value_in=root)
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    !! Handle base cases, or call to general implementation.
    IF (root .EQ. 1) THEN
       CALL Invert(InputMat, OutputMat, solver_parameters)
    ELSE IF (root .EQ. 2) THEN
       CALL InverseSquareRoot(InputMat, OutputMat, solver_parameters)
    ELSE IF (root .EQ. 3) THEN
       WRITE(*,*) "NOT IMPLEMENTED YET!"
    ELSE IF (root .EQ. 4) THEN
       CALL SquareRoot(InputMat, TempMat, solver_parameters)
       CALL InverseSquareRoot(TempMat, OutputMat, solver_parameters)
       CALL DestructDistributedSparseMatrix(TempMat)
    ELSE
       CALL ComputeInverseRootImplemention(InputMat, OutputMat, root, &
            & solver_parameters)
    END IF

    IF (solver_parameters%be_verbose .AND. IsRoot()) THEN
       CALL ExitSubLog
    END IF
  END SUBROUTINE ComputeInverseRoot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a general inverse matrix root for root > 4.
  !! Implementation details...
  SUBROUTINE ComputeInverseRootImplemention(InputMat, OutputMat, root, &
       & solver_parameters_in)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in)  :: InputMat
    TYPE(DistributedSparseMatrix), INTENT(inout) :: OutputMat
    INTEGER, INTENT(in) :: root
    TYPE(IterativeSolverParameters), INTENT(in), OPTIONAL :: solver_parameters_in
    REAL(NTREAL), PARAMETER :: NEGATIVE_ONE = -1.0
    !! Handling Solver Parameters
    TYPE(IterativeSolverParameters) :: solver_parameters
    !! Local Matrices
    TYPE(DistributedSparseMatrix) :: SqrtMat, FthrtMat
    TYPE(DistributedSparseMatrix) :: IdentityMat
    TYPE(DistributedSparseMatrix) :: Mk
    TYPE(DistributedSparseMatrix) :: IntermediateMat
    TYPE(DistributedSparseMatrix) :: IntermediateMatP
    TYPE(DistributedSparseMatrix) :: Temp
    !! Local Variables
    INTEGER :: target_root
    REAL(NTREAL) :: e_min
    REAL(NTREAL) :: e_max
    REAL(NTREAL) :: scaling_factor
    REAL(NTREAL) :: norm_value
    !! Temporary Variables
    INTEGER :: min_size, max_size
    REAL(NTREAL) :: sparsity
    INTEGER :: outer_counter
    INTEGER :: inner_counter
    TYPE(DistributedMatrixMemoryPool_t) :: pool

    !! Handle The Optional Parameters
    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters()
    END IF

    IF (solver_parameters_in%be_verbose .AND. IsRoot()) THEN
       CALL WriteHeader("Root Solver")
       CALL EnterSubLog
       CALL WriteCitation("nicholas2008functions")
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    !! Compute The Scaling Factor
    CALL EigenCircle(InputMat, e_min, e_max)
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
    CALL SquareRoot(InputMat, SqrtMat, solver_parameters)
    CALL SquareRoot(SqrtMat, FthrtMat, solver_parameters)
    CALL DestructDistributedSparseMatrix(SqrtMat)

    !!! Setup the Matrices
    CALL ConstructEmpty(IdentityMat, InputMat%actual_matrix_dimension)
    CALL FillDistributedIdentity(IdentityMat)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(FthrtMat, FthrtMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
       CALL PermuteMatrix(IdentityMat, IdentityMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF

    CALL CopyDistributedSparseMatrix(IdentityMat, OutputMat)
    CALL ScaleDistributedSparseMatrix(OutputMat, 1.0/scaling_factor)

    CALL CopyDistributedSparseMatrix(FthrtMat, Mk)
    CALL ScaleDistributedSparseMatrix(Mk, 1.0/(scaling_factor**target_root))
    CALL DestructDistributedSparseMatrix(FthrtMat)

    CALL ConstructEmpty(IntermediateMat, InputMat%actual_matrix_dimension)
    CALL ConstructEmpty(IntermediateMatP, InputMat%actual_matrix_dimension)
    CALL ConstructEmpty(Temp, InputMat%actual_matrix_dimension)

    outer_counter = 1
    IF (solver_parameters%be_verbose .AND. IsRoot()) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    norm_value = solver_parameters%converge_diff + 1.0d+0
    DO outer_counter = 1,solver_parameters%max_iterations
       IF (solver_parameters%be_verbose .AND. IsRoot() .AND. &
            & outer_counter .GT. 1) THEN
          CALL WriteListElement(key="Round", int_value_in=outer_counter-1)
          CALL EnterSubLog
          CALL WriteListElement(key="Convergence", float_value_in=norm_value)
          CALL ExitSubLog
       END IF

       CALL CopyDistributedSparseMatrix(IdentityMat, IntermediateMat)
       CALL ScaleDistributedSparseMatrix(IntermediateMat, &
            & REAL(target_root+1,NTREAL))
       CALL IncrementDistributedSparseMatrix(Mk, IntermediateMat, &
            & alpha_in=NEGATIVE_ONE)
       CALL ScaleDistributedSparseMatrix(IntermediateMat, &
            & REAL(1.0,NTREAL)/target_root)

       CALL DistributedGemm(OutputMat, IntermediateMat, Temp, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       CALL CopyDistributedSparseMatrix(Temp, OutputMat)

       CALL CopyDistributedSparseMatrix(IntermediateMat, IntermediateMatP)
       DO inner_counter = 1, target_root-1
          CALL DistributedGemm(IntermediateMat, IntermediateMatP, Temp, &
               & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
          CALL CopyDistributedSparseMatrix(Temp, IntermediateMatP)
       END DO

       CALL DistributedGemm(IntermediateMatP, Mk, Temp, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       CALL CopyDistributedSparseMatrix(Temp, Mk)

       CALL IncrementDistributedSparseMatrix(IdentityMat, Temp, &
            & alpha_in=NEGATIVE_ONE)
       norm_value = DistributedSparseNorm(Temp)

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    IF (solver_parameters%be_verbose) THEN
       CALL GetLoadBalance(OutputMat,min_size,max_size)
       sparsity = GetSize(OutputMat)/ &
            & (OutputMat%actual_matrix_dimension**2)
       IF (IsRoot()) THEN
          CALL ExitSubLog
          CALL WriteElement(key="Total_Iterations",int_value_in=outer_counter-1)
          CALL WriteHeader("Load_Balance")
          CALL EnterSubLog
          CALL WriteListElement(key="min_size", int_value_in=min_size)
          CALL WriteListElement(key="max_size", int_value_in=max_size)
          CALL ExitSubLog
          CALL WriteElement(key="Sparsity", float_value_in=sparsity)
       END IF
    END IF

    IF (MOD(root,4) .EQ. 1 .OR. MOD(root,4) .EQ. 3) THEN
       CALL DistributedGemm(OutputMat, OutputMat, Temp, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       CALL DistributedGemm(Temp, Temp, OutputMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    ELSE IF (MOD(root,4) .NE. 0) THEN
       CALL DistributedGemm(OutputMat, OutputMat, Temp, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       CALL CopyDistributedSparsematrix(Temp, OutputMat)
    END IF

    !! Undo Load Balancing Step
    CALL StartTimer("Load Balance")
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(OutputMat, OutputMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool)
    END IF
    CALL StopTimer("Load Balance")

    !!! Cleanup
    IF (solver_parameters_in%be_verbose .AND. IsRoot()) THEN
       CALL ExitSubLog
    END IF
    CALL DestructDistributedSparseMatrix(IdentityMat)
    CALL DestructDistributedSparseMatrix(Mk)
    CALL DestructDistributedSparseMatrix(IntermediateMat)
    CALL DestructDistributedSparseMatrix(IntermediateMatP)
    CALL DestructDistributedSparseMatrix(Temp)
  END SUBROUTINE ComputeInverseRootImplemention
END MODULE RootSolversModule
