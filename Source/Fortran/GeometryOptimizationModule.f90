!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Geometry Optimization
MODULE GeometryOptimizationModule
  USE DataTypesModule
  USE MatrixMemoryPoolPModule
  USE MatrixPSAlgebraModule
  USE MatrixPSModule
  USE IterativeSolversModule
  USE LoadBalancerModule
  USE LoggingModule
  USE SquareRootSolversModule
  USE TimerModule
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Solvers
  PUBLIC :: PurificationExtrapolate
  PUBLIC :: LowdinExtrapolate
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a new guess at the Density Matrix after updating the geometry.
  !! Based on the purification algorithm in \cite niklasson2010trace .
  !! @param[in] PreviousDensity to extrapolate from.
  !! @param[in] Overlap the overlap matrix of the new geometry.
  !! @param[in] nel the number of electrons.
  !! @param[out] NewDensity the extrapolated density.
  !! @param[in] solver_parameters_in parameters for the solver
  SUBROUTINE PurificationExtrapolate(PreviousDensity, Overlap, nel, NewDensity,&
       & solver_parameters_in)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN)  :: PreviousDensity, Overlap
    INTEGER, INTENT(IN) :: nel
    TYPE(Matrix_ps), INTENT(INOUT) :: NewDensity
    TYPE(IterativeSolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(Matrix_ps) :: WorkingDensity
    TYPE(Matrix_ps) :: WorkingOverlap
    TYPE(Matrix_ps) :: AddBranch, SubtractBranch
    TYPE(Matrix_ps) :: TempMat1, TempMat2
    TYPE(Matrix_ps) :: Identity
    !! Local Variables
    REAL(NTREAL), PARAMETER :: NEGATIVE_ONE = -1.0
    REAL(NTREAL) :: subtract_trace
    REAL(NTREAL) :: add_trace
    REAL(NTREAL) :: trace_value
    REAL(NTREAL) :: norm_value
    !! Temporary Variables
    TYPE(MatrixMemoryPool_p) :: pool1
    INTEGER :: outer_counter
    INTEGER :: total_iterations

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Density Matrix Extrapolator")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="Purification")
       CALL WriteCitation("niklasson2010trace")
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrixDS(NewDensity, &
         & PreviousDensity%actual_matrix_dimension, &
         & PreviousDensity%process_grid)
    CALL ConstructEmptyMatrixDS(WorkingDensity, &
         & PreviousDensity%actual_matrix_dimension, &
         & PreviousDensity%process_grid)
    CALL ConstructEmptyMatrixDS(WorkingOverlap, &
         & PreviousDensity%actual_matrix_dimension, &
         & PreviousDensity%process_grid)
    CALL ConstructEmptyMatrixDS(Identity, &
         & PreviousDensity%actual_matrix_dimension, &
         & PreviousDensity%process_grid)
    CALL FillDistributedIdentity(Identity)

    !! Compute the working hamiltonian.
    CALL CopyMatrixDS(PreviousDensity, WorkingDensity)
    CALL CopyMatrixDS(Overlap, WorkingOverlap)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(WorkingDensity, WorkingDensity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
       CALL PermuteMatrix(WorkingOverlap, WorkingOverlap, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
    END IF

    !! Finish Setup
    CALL CopyMatrixDS(WorkingDensity, NewDensity)
    CALL CopyMatrixDS(WorkingDensity, AddBranch)
    CALL CopyMatrixDS(WorkingDensity, SubtractBranch)

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0d+0
    DO outer_counter = 1,solver_parameters%max_iterations
       !! Figure Out Sigma Value. After which, XnS is stored in TempMat
       IF (outer_counter .GT. 1) THEN
          CALL DistributedGemm(AddBranch, WorkingOverlap, TempMat1, &
               & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
          add_trace = Trace(TempMat1)
          CALL DistributedGemm(SubtractBranch, WorkingOverlap, TempMat2, &
               & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
          subtract_trace = Trace(TempMat2)
          IF (ABS(nel - add_trace) .GT. ABS(nel - subtract_trace)) THEN
             !! Subtract Branch
             trace_value = subtract_trace
             CALL IncrementDistributedSparseMatrix(AddBranch, NewDensity, &
                  & NEGATIVE_ONE)
             norm_value = DistributedSparseNorm(NewDensity)
             CALL CopyMatrixDS(AddBranch, NewDensity)
             CALL CopyMatrixDS(TempMat2, TempMat1)
          ELSE
             !! Add Branch
             trace_value = add_trace
             CALL IncrementDistributedSparseMatrix(SubtractBranch, NewDensity, &
                  & NEGATIVE_ONE)
             norm_value = DistributedSparseNorm(NewDensity)
             CALL CopyMatrixDS(SubtractBranch, NewDensity)
          END IF
       ELSE
          CALL DistributedGemm(NewDensity, WorkingOverlap, TempMat1, &
               & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
       END IF

       IF (solver_parameters%be_verbose .AND. outer_counter .GT. 1) THEN
          CALL WriteListElement(key="Round", int_value_in=outer_counter-1)
          CALL EnterSubLog
          CALL WriteListElement(key="Convergence", float_value_in=norm_value)
          CALL WriteListElement(key="Trace", float_value_in=trace_value)
          CALL WriteListElement(key="AddTrace", float_value_in=add_trace)
          CALL WriteListElement(key="SubtractTrace", float_value_in=subtract_trace)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF

       !! Compute (I - XnS)Xn
       CALL IncrementDistributedSparseMatrix(Identity, TempMat1, &
            & alpha_in=NEGATIVE_ONE)
       CALL ScaleDistributedSparseMatrix(TempMat1, NEGATIVE_ONE)
       CALL DistributedGemm(TempMat1, NewDensity, TempMat2, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)

       !! Subtracted Version Xn - (I - XnS)Xn
       CALL CopyMatrixDS(NewDensity, SubtractBranch)
       CALL IncrementDistributedSparseMatrix(TempMat2, SubtractBranch, &
            & NEGATIVE_ONE)

       !! Added Version Xn + (I - XnS)Xn
       CALL CopyMatrixDS(NewDensity, AddBranch)
       CALL IncrementDistributedSparseMatrix(TempMat2, AddBranch)

    END DO
    total_iterations = outer_counter-1
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations",int_value_in=total_iterations)
       CALL PrintMatrixInformation(NewDensity)
    END IF

    !! Undo Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(NewDensity, NewDensity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
    END IF

    !! Cleanup
    CALL DestructMatrixDS(WorkingDensity)
    CALL DestructMatrixDS(WorkingOverlap)
    CALL DestructMatrixDS(Identity)
    CALL DestructMatrixDS(TempMat1)
    CALL DestructMatrixDS(TempMat2)
    CALL DestructMatrixDS(AddBranch)
    CALL DestructMatrixDS(SubtractBranch)
    CALL DestructMatrixMemoryPoolD(pool1)

    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

  END SUBROUTINE PurificationExtrapolate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a new guess at the Density Matrix after updating the geometry.
  !! Based on the lowdin algorithm in \cite exner2002comparison .
  !! @param[in] PreviousDensity to extrapolate from.
  !! @param[in] OldOverlap the overlap matrix of the old geometry.
  !! @param[in] NewOverlap the overlap matrix of the new geometry.
  !! @param[out] NewDensity the extrapolated density.
  !! @param[in] solver_parameters_in parameters for the solver
  SUBROUTINE LowdinExtrapolate(PreviousDensity, OldOverlap, NewOverlap, &
       & NewDensity, solver_parameters_in)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN)  :: PreviousDensity
    TYPE(Matrix_ps), INTENT(IN)  :: OldOverlap
    TYPE(Matrix_ps), INTENT(IN)  :: NewOverlap
    TYPE(Matrix_ps), INTENT(INOUT) :: NewDensity
    TYPE(IterativeSolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(IterativeSolverParameters_t) :: solver_parameters
    !! Local Matrices
    TYPE(Matrix_ps) :: SQRMat
    TYPE(Matrix_ps) :: ISQMat
    TYPE(Matrix_ps) :: TempMat
    !! Temporary Variables
    TYPE(MatrixMemoryPool_p) :: pool1

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = IterativeSolverParameters_t()
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Density Matrix Extrapolator")
       CALL EnterSubLog
       CALL WriteElement(key="Method", text_value_in="Lowdin")
       CALL WriteCitation("exner2002comparison")
       CALL PrintIterativeSolverParameters(solver_parameters)
    END IF

    CALL SquareRoot(OldOverlap, SQRMat, solver_parameters)
    CALL InverseSquareRoot(NewOverlap, ISQMat, solver_parameters)

    CALL DistributedGemm(SQRMat, PreviousDensity, TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
    CALL DistributedGemm(TempMat, SQRMat, NewDensity, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
    CALL DistributedGemm(ISQMat, NewDensity, TempMat, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
    CALL DistributedGemm(TempMat, ISQMat, NewDensity, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)

    CALL DestructMatrixDS(SQRMat)
    CALL DestructMatrixDS(ISQMat)
    CALL DestructMatrixDS(TempMat)

    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

  END SUBROUTINE LowdinExtrapolate
END MODULE GeometryOptimizationModule
