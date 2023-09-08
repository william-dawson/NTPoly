!> A module for monitoring convergence of an iterative algorithim.
!! In basic mode, we monitor that the last value isn't below the tight cutoff.
!! In automatic mode we monitor the following conditions:
!! o The last value can't be negative
!! o The moving average is within an order of magnitude
!! o The value isn't above the loose cutoff
MODULE ConvergenceMonitor
  USE DataTypesModule, ONLY : NTREAL
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteElement, &
       & WriteListElement
  IMPLICIT NONE
  PRIVATE
  !> Monitor convergence with a moving average
  TYPE, PUBLIC :: Monitor_t
     !> The small window
     REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: win_short
     !> The large window
     REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: win_long
     !> The number of values that have been added
     INTEGER :: nval
     !> We aren't converged if the average isn't below this.
     REAL(NTREAL) :: loose_cutoff
     !> We definitely are converged if the last value is below this.
     REAL(NTREAL) :: tight_cutoff
     !> Whether to do automatic exit
     LOGICAL :: automatic
  END TYPE Monitor_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructMonitor
  PUBLIC :: DestructMonitor
  PUBLIC :: AppendValue
  PUBLIC :: CheckConverged
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct the monitor.
  SUBROUTINE ConstructMonitor(this, short_len_in, long_len_in, &
       & loose_cutoff_in, tight_cutoff_in, automatic_in)
    !> The monitor to construct.
    TYPE(Monitor_t), INTENT(INOUT) :: this
    !> The length of the short window (default: 3)
    INTEGER, INTENT(IN), OPTIONAL :: short_len_in
    !> The length of the long window (default: 6)
    INTEGER, INTENT(IN), OPTIONAL :: long_len_in
    !> If the average is greater than this than we aren't 
    !! converged (default: 0.001)
    REAL(NTREAL), INTENT(IN), OPTIONAL :: loose_cutoff_in
    !> If the last value is less than this, we definitely are converged 
    !! (default: 1e-8).
    REAL(NTREAL), INTENT(IN), OPTIONAL :: tight_cutoff_in
    !> Whether to use automatic mode (default: .true.)
    LOGICAL, INTENT(IN), OPTIONAL :: automatic_in

    CALL DestructMonitor(this)

    !> Allocate Windows
    IF (PRESENT(short_len_in)) THEN
       ALLOCATE(this%win_short(short_len_in))
    ELSE
       ALLOCATE(this%win_short(3))
    END IF
    this%win_short = 0

    IF (PRESENT(long_len_in)) THEN
       ALLOCATE(this%win_long(long_len_in))
    ELSE
       ALLOCATE(this%win_long(6))
    END IF
    this%win_long = 0

    IF (PRESENT(loose_cutoff_in)) THEN
       this%loose_cutoff = loose_cutoff_in
    ELSE
       this%loose_cutoff = 1E-3_NTREAL
    END IF

    IF (PRESENT(tight_cutoff_in)) THEN
       this%tight_cutoff = tight_cutoff_in
    ELSE
       this%tight_cutoff = 1E-8_NTREAL
    END IF

    IF (PRESENT(automatic_in)) THEN
       this%automatic = automatic_in
    ELSE
       this%automatic = .TRUE.
    END IF

    this%nval = 0
  END SUBROUTINE ConstructMonitor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct the monitor.
  SUBROUTINE DestructMonitor(this)
    !> The monitor to destruct.
    TYPE(Monitor_t), INTENT(INOUT) :: this
    IF (ALLOCATED(this%win_short)) DEALLOCATE(this%win_short)
    IF (ALLOCATED(this%win_long)) DEALLOCATE(this%win_long)
  END SUBROUTINE DestructMonitor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Add a value to the window
  SUBROUTINE AppendValue(this, val)
    !> Add a value to the window
    TYPE(Monitor_t), INTENT(INOUT) :: this
    !> Value to add
    REAL(NTREAL), INTENT(IN) :: val
    !! Local Variables
    INTEGER :: II

    !! Shift
    DO II = SIZE(this%win_short), 2, -1
       this%win_short(II - 1) = this%win_short(II) 
    END DO
    DO II = SIZE(this%win_long), 2, -1
       this%win_long(II - 1) = this%win_long(II)
    END DO

    !! Append
    this%win_short(SIZE(this%win_short)) = val
    this%win_long(SIZE(this%win_long)) = val
    this%nval = this%nval + 1
  END SUBROUTINE AppendValue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Check if Convergence Has Been Achieved
  FUNCTION CheckConverged(this, be_verbose) RESULT(conv)
    !> The window to test
    TYPE(Monitor_t), INTENT(IN) :: this
    !> True if converged.
    LOGICAL :: conv
    !> True if we should print convergence info
    LOGICAL :: be_verbose
    !! Local Variables
    REAL(NTREAL) :: avg_short, avg_long
    REAL(NTREAL) :: last, last2

    !! Basic Check
    last = this%win_short(SIZE(this%win_short))
    last2 = this%win_short(SIZE(this%win_short) - 1)
    IF (be_verbose) THEN
       CALL WriteListElement(key = "Convergence", VALUE = last)
    END IF
    IF (ABS(last) .GT. this%tight_cutoff) THEN
       conv = .FALSE.
    ELSE
       conv = .TRUE.
       CALL EnterSubLog
       CALL WriteElement(key = "Trigger", VALUE = "Tight Criteria")
       CALL ExitSubLog
    END IF
    IF (.NOT. this%automatic .OR. conv) RETURN

    !! Automatic disabled
    conv = .TRUE.

    !! First check that we've seen enough values to make a judgement.
    IF (this%nval .LT. SIZE(this%win_long)) conv = .FALSE.

    !! Now check that the two windows are within an order of magnitude
    avg_short = SUM(this%win_short) / SIZE(this%win_short)
    avg_long = SUM(this%win_long) / SIZE(this%win_long)
    IF (.NOT. (10 * avg_short .GT. avg_long .AND. &
         & avg_short / 10 .LT. avg_long)) THEN
       conv = .FALSE.
    END IF
    IF (be_verbose) THEN
       CALL EnterSubLog
       CALL WriteElement(key = "Avg Short", VALUE = avg_short)
       CALL WriteElement(key = "Avg Long", VALUE = avg_long)
       CALL ExitSubLog
    END IF

    !! No convergence if the lastest value is negative
    IF (last .LT. 0) conv = .FALSE.

    !! No convergence if the absolute value is decreasing in magnitude
    IF (ABS(last) < ABS(last2)) conv = .FALSE.

    !! No convergence if the average value is greater than the cutoff
    IF (avg_long .GT. this%loose_cutoff) conv = .FALSE.

    IF (conv) THEN
       CALL EnterSubLog
       CALL WriteElement(key = "Trigger", VALUE = "Automatic")
       CALL ExitSubLog
    END IF
  END FUNCTION CheckConverged
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE ConvergenceMonitor
