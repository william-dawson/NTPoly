!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Storing The Parameters For Iterative Solvers.
MODULE SolverParametersModule
  USE ConvergenceMonitor, ONLY : Monitor_t, DestructMonitor
  USE DataTypesModule, ONLY : NTREAL
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteElement, &
       & WriteHeader
  USE PermutationModule, ONLY : Permutation_t, CopyPermutation, &
       & DestructPermutation
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A class for passing parameters to an iterative solver.
  TYPE, PUBLIC :: SolverParameters_t
     !> When do we consider a calculation converged.
     REAL(NTREAL) :: converge_diff
     !> Maximum number of iterations of a solver before termination.
     INTEGER :: max_iterations
     !> Threshold for sparse multiplication and addition.
     REAL(NTREAL) :: threshold
     !> If true, the sparse solver prints out information each loop iteration.
     LOGICAL :: be_verbose
     !> If true, the sparse solver will try and load balance before calculation.
     LOGICAL :: do_load_balancing
     !> The permutation used for load balancing.
     TYPE(Permutation_t) :: BalancePermutation
     !> Thresholds for step size searches.
     REAL(NTREAL) :: step_thresh
     !> Whether to do an automatic convergence detection.
     LOGICAL :: monitor_convergence
     !> The convergence monitor
     TYPE(Monitor_t) :: monitor
  END TYPE SolverParameters_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructSolverParameters
  PUBLIC :: CopySolverParameters
  PUBLIC :: SetParametersConvergeDiff
  PUBLIC :: SetParametersMaxIterations
  PUBLIC :: SetParametersThreshold
  PUBLIC :: SetParametersBeVerbose
  PUBLIC :: SetParametersLoadBalance
  PUBLIC :: SetParametersStepThreshold
  PUBLIC :: SetParametersMonitorConvergence
  PUBLIC :: PrintParameters
  PUBLIC :: DestructSolverParameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The default convergence difference.
  REAL(NTREAL), PARAMETER, PUBLIC :: CONVERGENCE_DIFF_CONST = 1e-6_NTREAL
  !> The default maximum number of iterations.
  INTEGER, PARAMETER, PUBLIC :: MAX_ITERATIONS_CONST = 1000
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a data type which stores iterative solver parameters.
  SUBROUTINE ConstructSolverParameters(this, converge_diff_in, threshold_in, &
       & max_iterations_in, be_verbose_in, BalancePermutation_in, &
       & step_thresh_in, monitor_convergence_in)
    !> The parameters to construct.
    TYPE(SolverParameters_t), INTENT(INOUT) :: this
    !> Converge_diff_in the difference between iterations to consider
    !> a calculation converged.
    REAL(NTREAL), INTENT(IN), OPTIONAL :: converge_diff_in
    !> The zero threshold
    REAL(NTREAL), INTENT(IN), OPTIONAL :: threshold_in
    !> The maximum number of iterations to perform
    INTEGER, INTENT(IN), OPTIONAL :: max_iterations_in
    !> Whether to print during the calculation (default = False)
    LOGICAL, INTENT(IN), OPTIONAL :: be_verbose_in
    !> For load balancing
    TYPE(Permutation_t), INTENT(IN), OPTIONAL :: BalancePermutation_in
    !> Step size for differential equation solvers.
    REAL(NTREAL), INTENT(IN), OPTIONAL :: step_thresh_in
    !> Whether to do automatic convergence monitoring (default = True).
    LOGICAL, INTENT(IN), OPTIONAL :: monitor_convergence_in

    CALL DestructSolverParameters(this)

    !! Optional Parameters
    IF (.NOT. PRESENT(converge_diff_in)) THEN
       this%converge_diff = CONVERGENCE_DIFF_CONST
    ELSE
       this%converge_diff = converge_diff_in
    END IF
    IF (.NOT. PRESENT(threshold_in)) THEN
       this%threshold = 0.0_NTREAL
    ELSE
       this%threshold = threshold_in
    END IF
    IF (.NOT. PRESENT(max_iterations_in)) THEN
       this%max_iterations = MAX_ITERATIONS_CONST
    ELSE
       this%max_iterations = max_iterations_in
    END IF
    IF (.NOT. PRESENT(be_verbose_in)) THEN
       this%be_verbose = .FALSE.
    ELSE
       this%be_verbose = be_verbose_in
    END IF
    IF (.NOT. PRESENT(BalancePermutation_in)) THEN
       this%do_load_balancing = .FALSE.
    ELSE
       this%do_load_balancing = .TRUE.
       CALL CopyPermutation(BalancePermutation_in, this%BalancePermutation)
    END IF
    IF (.NOT. PRESENT(step_thresh_in)) THEN
       this%step_thresh = 1E-2_NTREAL
    ELSE
       this%step_thresh = step_thresh_in
    END IF
    IF (.NOT. PRESENT(monitor_convergence_in)) THEN
       this%monitor_convergence = .TRUE.
    ELSE
       this%monitor_convergence = monitor_convergence_in
    END IF
  END SUBROUTINE ConstructSolverParameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE CopySolverParameters(paramA, paramB)
    !> Parameters to copy
    TYPE(SolverParameters_t), INTENT(IN) :: paramA
    !> paramB = paramA
    TYPE(SolverParameters_t), INTENT(INOUT) :: paramB

    CALL DestructSolverParameters(paramB)
    paramB = paramA
  END SUBROUTINE CopySolverParameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the convergence difference.
  PURE SUBROUTINE SetParametersConvergeDiff(this, new_value)
    !> The parameter object.
    TYPE(SolverParameters_t), INTENT(INOUT) :: this
    !> Value to set it to.
    REAL(NTREAL), INTENT(IN) :: new_value

    this%converge_diff = new_value
  END SUBROUTINE SetParametersConvergeDiff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the max iterations.
  PURE SUBROUTINE SetParametersMaxIterations(this, new_value)
    !> The parameter object.
    TYPE(SolverParameters_t), INTENT(INOUT) :: this
    !> Value to set it to.
    INTEGER, INTENT(IN) :: new_value

    this%max_iterations = new_value
  END SUBROUTINE SetParametersMaxIterations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the threshold.
  PURE SUBROUTINE SetParametersThreshold(this, new_value)
    !> The parameter object.
    TYPE(SolverParameters_t), INTENT(INOUT) :: this
    !> Value to set it to.
    REAL(NTREAL), INTENT(IN) :: new_value

    this%threshold = new_value
  END SUBROUTINE SetParametersThreshold
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the verbosity.
  PURE SUBROUTINE SetParametersBeVerbose(this, new_value)
    !> The parameter object.
    TYPE(SolverParameters_t), INTENT(INOUT) :: this
    !> Value to set it to.
    LOGICAL, INTENT(IN) :: new_value

    this%be_verbose = new_value
  END SUBROUTINE SetParametersBeVerbose
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the load balance.
  PURE SUBROUTINE SetParametersLoadBalance(this, new_value)
    !> The parameter object.
    TYPE(SolverParameters_t), INTENT(INOUT) :: this
    !> Value to set it to.
    TYPE(Permutation_t), INTENT(IN) :: new_value

    this%do_load_balancing = .TRUE.
    this%BalancePermutation = new_value
  END SUBROUTINE SetParametersLoadBalance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the step threshold.
  PURE SUBROUTINE SetParametersStepThreshold(this, new_value)
    !> The parameter object.
    TYPE(SolverParameters_t), INTENT(INOUT) :: this
    !> Value to set it to.
    REAL(NTREAL), INTENT(IN) :: new_value

    this%step_thresh = new_value
  END SUBROUTINE SetParametersStepThreshold
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the monitoring of convergence
  PURE SUBROUTINE SetParametersMonitorConvergence(this, new_value)
    !> The parameter object.
    TYPE(SolverParameters_t), INTENT(INOUT) :: this
    !> Value to set it to.
    LOGICAL, INTENT(IN) :: new_value

    this%monitor_convergence = new_value
  END SUBROUTINE SetParametersMonitorConvergence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out the iterative solver parameter values.
  SUBROUTINE PrintParameters(this)
    !> The parameter object.
    TYPE(SolverParameters_t), INTENT(IN) :: this

    CALL WriteHeader("Solver Parameters")
    CALL EnterSubLog
    CALL WriteElement(key = "Verbosity", &
         & VALUE = this%be_verbose)
    CALL WriteElement(key = "Load Balancing", &
         & VALUE = this%do_load_balancing)
    CALL WriteElement(key = "Convergence Difference", &
         & VALUE = this%converge_diff)
    CALL WriteElement(key = "Threshold", &
         & VALUE = this%threshold)
    CALL WriteElement(key = "Maximum Iterations", &
         & VALUE = this%max_iterations)
    CALL WriteElement(key = "Step Threshold", &
         & VALUE = this%step_thresh)
    CALL WriteElement(key = "Monitor Convergence", &
         & VALUE = this%monitor_convergence)
    CALL ExitSubLog
  END SUBROUTINE PrintParameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Cleanup the solver parameters datastructure.
  PURE SUBROUTINE DestructSolverParameters(this)
    !> The parameter object.
    TYPE(SolverParameters_t), INTENT(INOUT) :: this

    CALL DestructPermutation(this%BalancePermutation)
    CALL DestructMonitor(this%monitor)
  END SUBROUTINE DestructSolverParameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE SolverParametersModule
