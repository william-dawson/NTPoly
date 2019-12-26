!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module to do timings.
MODULE TimerModule
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteElement, &
       & WriteHeader
  USE ProcessGridModule, ONLY : global_grid
  USE NTMPIModule
  IMPLICIT NONE
  PRIVATE
  LOGICAL :: is_initialized = .FALSE.
  CHARACTER(len=20), DIMENSION(:), ALLOCATABLE :: timer_list
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: start_times
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: elapsed_times
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: RegisterTimer
  PUBLIC :: StartTimer
  PUBLIC :: StopTimer
  PUBLIC :: PrintTimer
  PUBLIC :: PrintAllTimers
  PUBLIC :: PrintAllTimersDistributed
CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Register a timer with the timer module.  Call this before using that timer.
  SUBROUTINE RegisterTimer(timer_name)
    !> Name of the timer.
    CHARACTER(len=*), INTENT(IN) :: timer_name
    !! Local Data
    CHARACTER(len=20), DIMENSION(:), ALLOCATABLE :: temp_timer_list
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: temp_start_times
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: temp_elapsed_times

    IF (is_initialized) THEN
       ALLOCATE(temp_timer_list(SIZE(timer_list)+1))
       ALLOCATE(temp_start_times(SIZE(start_times)+1))
       ALLOCATE(temp_elapsed_times(SIZE(elapsed_times)+1))
       temp_timer_list(:SIZE(timer_list)) = timer_list
       temp_start_times(:SIZE(start_times)) = start_times
       temp_elapsed_times(:SIZE(elapsed_times)) = elapsed_times
       CALL MOVE_ALLOC(temp_timer_list,timer_list)
       CALL MOVE_ALLOC(temp_start_times,start_times)
       CALL MOVE_ALLOC(temp_elapsed_times,elapsed_times)
       timer_list(SIZE(timer_list)) = timer_name
       elapsed_times(SIZE(timer_list)) = 0
    ELSE
       ALLOCATE(timer_list(1))
       ALLOCATE(start_times(1))
       ALLOCATE(elapsed_times(1))
       timer_list(1) = timer_name
       elapsed_times(1) = 0
       is_initialized = .TRUE.
    END IF
  END SUBROUTINE RegisterTimer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Start the clock running for a given timer.
  SUBROUTINE StartTimer(timer_name)
    !> Name of the timer. Must be registered.
    CHARACTER(len=*), INTENT(IN) :: timer_name
    !! Local Data
    INTEGER :: timer_position
    DOUBLE PRECISION :: temp_time

    temp_time = MPI_WTIME()
    timer_position = GetTimerPosition(timer_name)
    IF (timer_position > 0) THEN
       start_times(timer_position) = temp_time
    END IF
  END SUBROUTINE StartTimer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Stop the clock for a given timer.
  SUBROUTINE StopTimer(timer_name)
    !> Name of the timer. Must be registered.
    CHARACTER(len=*), INTENT(IN) :: timer_name
    !! Local Data
    INTEGER :: timer_position
    DOUBLE PRECISION :: temp_elapsed_time
    DOUBLE PRECISION :: temp_start_time

    timer_position = GetTimerPosition(timer_name)
    IF (timer_position > 0) THEN
       temp_elapsed_time = MPI_WTIME()
       temp_start_time = start_times(timer_position)
       elapsed_times(timer_position) = elapsed_times(timer_position) + &
            & temp_elapsed_time - temp_start_time
    END IF
  END SUBROUTINE StopTimer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out the elapsed time for a given timer.
  SUBROUTINE PrintTimer(timer_name)
    !> Name of the timer. Must be registered.
    CHARACTER(len=*), INTENT(IN) :: timer_name
    !! Local Data
    INTEGER :: timer_position

    timer_position = GetTimerPosition(timer_name)
    CALL WriteHeader("Timers")
    CALL EnterSubLog
    IF (timer_position > 0) THEN
       CALL WriteElement(key=timer_name, &
            & value=elapsed_times(timer_position))
    END IF
    CALL ExitSubLog
  END SUBROUTINE PrintTimer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out the elapsed time for each timer on this process.
  SUBROUTINE PrintAllTimers()
    !! Local Data
    INTEGER :: timer_position

    CALL WriteHeader("Timers")
    CALL EnterSubLog
    DO timer_position = LBOUND(timer_list,dim=1), UBOUND(timer_list,dim=1)
       CALL WriteElement(key=timer_list(timer_position), &
            & value=elapsed_times(timer_position))
    END DO
    CALL ExitSubLog
  END SUBROUTINE PrintAllTimers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out the elapsed time for each timer based on the max value across
  !> processes.
  SUBROUTINE PrintAllTimersDistributed()
    !! Local Data
    INTEGER          :: timer_position
    DOUBLE PRECISION :: elapsed
    DOUBLE PRECISION :: max_time
    INTEGER          :: ierr

    CALL WriteHeader("Timers")
    CALL EnterSubLog

    DO timer_position = LBOUND(timer_list,dim=1), UBOUND(timer_list,dim=1)
       elapsed = elapsed_times(timer_position)
       CALL MPI_Allreduce(elapsed, max_time, 1, MPI_DOUBLE ,MPI_MAX, &
            & global_grid%global_comm, ierr)
       CALL WriteElement(key=timer_list(timer_position), &
            & value=max_time)
    END DO

    CALL ExitSubLog
  END SUBROUTINE PrintAllTimersDistributed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Figure out the position in the timer list where timer_name is.
  !> This is a utility routine.
  FUNCTION GetTimerPosition(timer_name) RESULT(timer_position)
    !! Parameters
    !> Name of the timer.
    CHARACTER(len=*), INTENT(IN) :: timer_name
    !> The position of the timer. 0 means the timer has not been registered.
    INTEGER :: timer_position
    !! Local Data
    INTEGER :: counter
    LOGICAL :: not_found

    not_found = .TRUE.

    IF (is_initialized) THEN
       DO counter=1, SIZE(timer_list)
          IF (timer_name .EQ. timer_list(counter)) THEN
             not_found = .FALSE.
             EXIT
          END IF
       END DO
    END IF

    IF (not_found) THEN
       timer_position = 0
    ELSE
       timer_position = counter
    END IF
  END FUNCTION GetTimerPosition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE TimerModule
