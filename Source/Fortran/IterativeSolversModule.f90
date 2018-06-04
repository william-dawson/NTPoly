!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Storing The Parameters For Iterative Solvers.
MODULE IterativeSolversModule
  USE DataTypesModule, ONLY : NTREAL
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteListElement, &
       & WriteHeader
  USE PermutationModule, ONLY : Permutation_t
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A class for passing parameters to an iterative solver.
  TYPE, PUBLIC :: IterativeSolverParameters_t
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
  END TYPE IterativeSolverParameters_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE IterativeSolverParameters_t
     MODULE PROCEDURE IterativeSolverParameters_init
  END INTERFACE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SetIterativeConvergeDiff
  PUBLIC :: SetIterativeMaxIterations
  PUBLIC :: SetIterativeThreshold
  PUBLIC :: SetIterativeBeVerbose
  PUBLIC :: SetIterativeLoadBalance
  PUBLIC :: PrintIterativeSolverParameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The default convergence difference.
  REAL(NTREAL), PARAMETER, PUBLIC :: CONVERGENCE_DIFF_CONST = 1e-6
  !> The default maximum number of iterations.
  INTEGER, PARAMETER, PUBLIC :: MAX_ITERATIONS_CONST = 100
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Handle the parameters
  !> Construct a data type which stores iterative solver parameters.
  !! @param[in] converge_diff_in the difference between iterations to consider
  !! a calculation converged (optional).
  !! @param[in] threshold_in the zero threshold (optional).
  !! @param[in] max_iterations_in the maximum number of iterations to perform
  !! (optional).
  !! @param[in] be_verbose_in whether to print during the calculation
  !! (optional, default = False).
  !! @param[in] BalancePermutation_in for load balancing (optional).
  PURE FUNCTION IterativeSolverParameters_init(converge_diff_in, threshold_in, &
       & max_iterations_in, be_verbose_in, BalancePermutation_in) RESULT(this)
    !! Parameters
    REAL(NTREAL), INTENT(IN), OPTIONAL :: converge_diff_in
    REAL(NTREAL), INTENT(IN), OPTIONAL :: threshold_in
    INTEGER, INTENT(IN), OPTIONAL :: max_iterations_in
    LOGICAL, INTENT(IN), OPTIONAL :: be_verbose_in
    TYPE(Permutation_t), INTENT(IN), OPTIONAL :: BalancePermutation_in
    TYPE(IterativeSolverParameters_t) :: this

    !! Optional Parameters
    IF (.NOT. PRESENT(converge_diff_in)) THEN
       this%converge_diff = CONVERGENCE_DIFF_CONST
    ELSE
       this%converge_diff = converge_diff_in
    END IF
    IF (.NOT. PRESENT(threshold_in)) THEN
       this%threshold = 0.0
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
       this%BalancePermutation = BalancePermutation_in
    END IF
  END FUNCTION IterativeSolverParameters_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the convergence difference.
  !! @param[inout] this the parameter object.
  !! @param[in] new_value to set it to.
  PURE SUBROUTINE SetIterativeConvergeDiff(this,new_value)
    TYPE(IterativeSolverParameters_t), INTENT(INOUT) :: this
    REAL(NTREAL), INTENT(IN) :: new_value

    this%converge_diff = new_value
  END SUBROUTINE SetIterativeConvergeDiff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the max iterations.
  !! @param[inout] this the parameter object.
  !! @param[in] new_value to set it to.
  PURE SUBROUTINE SetIterativeMaxIterations(this,new_value)
    TYPE(IterativeSolverParameters_t), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: new_value

    this%max_iterations = new_value
  END SUBROUTINE SetIterativeMaxIterations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the threshold.
  !! @param[inout] this the parameter object.
  !! @param[in] new_value to set it to.
  PURE SUBROUTINE SetIterativeThreshold(this,new_value)
    TYPE(IterativeSolverParameters_t), INTENT(INOUT) :: this
    REAL(NTREAL), INTENT(IN) :: new_value

    this%threshold = new_value
  END SUBROUTINE SetIterativeThreshold
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the verbosity.
  !! @param[inout] this the parameter object.
  !! @param[in] new_value to set it to.
  PURE SUBROUTINE SetIterativeBeVerbose(this,new_value)
    TYPE(IterativeSolverParameters_t), INTENT(INOUT) :: this
    LOGICAL, INTENT(IN) :: new_value

    this%be_verbose = new_value
  END SUBROUTINE SetIterativeBeVerbose
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the load balance.
  !! @param[inout] this the parameter object.
  !! @param[in] new_value to set it to.
  PURE SUBROUTINE SetIterativeLoadBalance(this,new_value)
    TYPE(IterativeSolverParameters_t), INTENT(INOUT) :: this
    TYPE(Permutation_t), INTENT(IN) :: new_value

    this%do_load_balancing = .TRUE.
    this%BalancePermutation = new_value
  END SUBROUTINE SetIterativeLoadBalance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out the iterative solver parameter values.
  !! @param[inout] this the parameter object.
  SUBROUTINE PrintIterativeSolverParameters(this)
    TYPE(IterativeSolverParameters_t), INTENT(IN) :: this

    CALL WriteHeader("Iterative Solver Parameters")
    CALL EnterSubLog
    CALL WriteListElement(key="be_verbose",bool_value_in=this%be_verbose)
    CALL WriteListElement(key="do_load_balancing", &
         & bool_value_in=this%do_load_balancing)
    CALL WriteListElement(key="converge_diff", &
         & float_value_in=this%converge_diff)
    CALL WriteListElement(key="threshold", &
         & float_value_in=this%threshold)
    CALL WriteListElement(key="max_iterations", &
         & int_value_in=this%max_iterations)
    CALL ExitSubLog
  END SUBROUTINE PrintIterativeSolverParameters
END MODULE IterativeSolversModule
