!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Parameters for a fixed size polynomial solver.
MODULE FixedSolversModule
  USE DataTypesModule, ONLY : NTREAL
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteListElement, &
       & WriteHeader
  USE PermutationModule, ONLY : Permutation_t
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A class for passing parameters to an iterative solver.
  TYPE, PUBLIC :: FixedSolverParameters_t
     !> Threshold for sparse multiplication and addition.
     REAL(NTREAL) :: threshold
     !> If true, the sparse solver prints out information each loop iteration.
     LOGICAL :: be_verbose
     !> If true, the sparse solver will try and load balance before calculation.
     LOGICAL :: do_load_balancing
     !> The permutation used for load balancing.
     TYPE(Permutation_t) :: BalancePermutation
  END TYPE FixedSolverParameters_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE FixedSolverParameters_t
     MODULE PROCEDURE FixedSolverParameters_init
  END INTERFACE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SetFixedThreshold
  PUBLIC :: SetFixedBeVerbose
  PUBLIC :: SetFixedLoadBalance
  PUBLIC :: PrintFixedSolverParameters
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a fixed solver parameteres data structure.
  !! @param[in] threshold_in the zero threshold value (optional, default=0).
  !! @param[in] be_verbose_in whether to print during the solve
  !! (optional, default=False)
  !! @param[in] BalancePermutation_in a permutation for load balancing.
  PURE FUNCTION FixedSolverParameters_init(threshold_in, &
       & be_verbose_in, BalancePermutation_in) RESULT(this)
    !! Parameters
    REAL(NTREAL), INTENT(IN), OPTIONAL :: threshold_in
    LOGICAL, INTENT(IN), OPTIONAL :: be_verbose_in
    TYPE(Permutation_t), INTENT(IN), OPTIONAL :: BalancePermutation_in
    TYPE(FixedSolverParameters_t) :: this

    !! Optional Parameters
    IF (.NOT. PRESENT(threshold_in)) THEN
       this%threshold = 0.0
    ELSE
       this%threshold = threshold_in
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
  END FUNCTION FixedSolverParameters_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the threshold.
  !! @param[inout] this the parameter object.
  !! @param[in] new_value to set it to.
  PURE SUBROUTINE SetFixedThreshold(this,new_value)
    TYPE(FixedSolverParameters_t), INTENT(INOUT) :: this
    REAL(NTREAL), INTENT(IN) :: new_value

    this%threshold = new_value
  END SUBROUTINE SetFixedThreshold
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the verbosity.
  !! @param[inout] this the parameter object.
  !! @param[in] new_value to set it to.
  PURE SUBROUTINE SetFixedBeVerbose(this,new_value)
    TYPE(FixedSolverParameters_t), INTENT(INOUT) :: this
    LOGICAL, INTENT(IN) :: new_value

    this%be_verbose = new_value
  END SUBROUTINE SetFixedBeVerbose
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the load balance.
  !! @param[inout] this the parameter object.
  !! @param[in] new_value to set it to.
  PURE SUBROUTINE SetFixedLoadBalance(this,new_value)
    TYPE(FixedSolverParameters_t), INTENT(INOUT) :: this
    TYPE(Permutation_t), INTENT(IN) :: new_value

    this%do_load_balancing = .TRUE.
    this%BalancePermutation = new_value
  END SUBROUTINE SetFixedLoadBalance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out the fixed solver parameter values.
  !! @param[inout] this the parameter object.
  SUBROUTINE PrintFixedSolverParameters(this)
    TYPE(FixedSolverParameters_t), INTENT(IN) :: this
    CALL WriteHeader("Fixed Solver Parameters")
    CALL EnterSubLog
    CALL WriteListElement(key="be_verbose",bool_value_in=this%be_verbose)
    CALL WriteListElement(key="threshold",float_value_in=this%threshold)
    CALL WriteListElement(key="do_load_balancing", &
         & bool_value_in=this%do_load_balancing)
    CALL ExitSubLog
  END SUBROUTINE PrintFixedSolverParameters
END MODULE FixedSolversModule
