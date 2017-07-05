!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Parameters for a fixed size polynomial solver.
MODULE FixedSolversModule
  USE DataTypesModule
  USE LoggingModule
  USE PermutationModule
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A class for passing parameters to an iterative solver.
  TYPE, PUBLIC :: FixedSolverParameters
     !> Threshold for sparse multiplication and addition.
     REAL(NTREAL) :: threshold
     !> If true, the sparse solver prints out information each loop iteration.
     LOGICAL :: be_verbose
     !> If true, the sparse solver will try and load balance before calculation.
     LOGICAL :: do_load_balancing
     !> The permutation used for load balancing.
     TYPE(Permutation_t) :: BalancePermutation
  END TYPE FixedSolverParameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE FixedSolverParameters
     MODULE PROCEDURE FixedSolverParameters_init
  END INTERFACE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SetFixedThreshold
  PUBLIC :: SetFixedBeVerbose
  PUBLIC :: SetFixedLoadBalance
  PUBLIC :: PrintFixedSolverParameters
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Handle the parameters
  PURE FUNCTION FixedSolverParameters_init(threshold_in, &
       & be_verbose_in, BalancePermutation_in) RESULT(this)
    !! Parameters
    REAL(NTREAL), INTENT(in), OPTIONAL :: threshold_in
    LOGICAL, INTENT(in), OPTIONAL :: be_verbose_in
    TYPE(Permutation_t), INTENT(in), OPTIONAL :: BalancePermutation_in
    TYPE(FixedSolverParameters) :: this

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
    TYPE(FixedSolverParameters), INTENT(inout) :: this
    REAL(NTREAL), INTENT(in) :: new_value

    this%threshold = new_value
  END SUBROUTINE SetFixedThreshold
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the verbosity.
  !! @param[inout] this the parameter object.
  !! @param[in] new_value to set it to.
  PURE SUBROUTINE SetFixedBeVerbose(this,new_value)
    TYPE(FixedSolverParameters), INTENT(inout) :: this
    LOGICAL, INTENT(in) :: new_value

    this%be_verbose = new_value
  END SUBROUTINE SetFixedBeVerbose
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of the load balance.
  !! @param[inout] this the parameter object.
  !! @param[in] new_value to set it to.
  PURE SUBROUTINE SetFixedLoadBalance(this,new_value)
    TYPE(FixedSolverParameters), INTENT(inout) :: this
    TYPE(Permutation_t), INTENT(in) :: new_value

    this%do_load_balancing = .TRUE.
    this%BalancePermutation = new_value
  END SUBROUTINE SetFixedLoadBalance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out the convergence values.
  !! @param[inout] this the parameter object.
  SUBROUTINE PrintFixedSolverParameters(this)
    TYPE(FixedSolverParameters), INTENT(in) :: this
    CALL WriteHeader("Fixed Solver Parameters")
    CALL EnterSubLog
    CALL WriteListElement(key="be_verbose",bool_value_in=this%be_verbose)
    CALL WriteListElement(key="threshold",float_value_in=this%threshold)
    CALL WriteListElement(key="do_load_balancing", &
         & bool_value_in=this%do_load_balancing)
    CALL ExitSubLog
  END SUBROUTINE PrintFixedSolverParameters
END MODULE FixedSolversModule
