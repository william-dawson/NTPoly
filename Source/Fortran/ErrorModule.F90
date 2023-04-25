!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module to do handle error passing.
MODULE ErrorModule
  USE NTMPIModule
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructError
  PUBLIC :: SetGenericError
  PUBLIC :: SetCustomError
  PUBLIC :: CheckMPIError
  PUBLIC :: CheckAllocError
  PUBLIC :: ErrorOccurred
  PUBLIC :: PrintError
  PUBLIC :: Cleanup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A type that can be passed around to accumulate errors.
  TYPE, PUBLIC :: Error_t
     PRIVATE
     !> Flag for whether or not an error has occurred.
     LOGICAL :: error_set
     !> Detailed description of the error.
     CHARACTER(len=1000) :: error_description
     !> Store a failed MPI call error.
     INTEGER :: mpi_error
     !> Flag for whether mpi error occurred.
     LOGICAL :: mpi_error_set
     !> Store a bad allocation call error.
     INTEGER :: alloc_error
     !> Flag for whether alloc error occurred.
     LOGICAL :: alloc_error_set
     !> Store a custom error.
     INTEGER :: custom_error
     !> Flag for whether a custom error occurred.
     LOGICAL :: custom_error_set
  END TYPE Error_t
CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Default constructor for an error type.
  SUBROUTINE ConstructError(this)
    !> The newly constructed error type
    TYPE(Error_t), INTENT(INOUT) :: this

    this%error_set = .FALSE.
    this%mpi_error_set = .FALSE.
    this%alloc_error_set = .FALSE.
  END SUBROUTINE ConstructError
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Routine to call if a generic error has occurred.
  SUBROUTINE SetGenericError(this, error_description, immediate_cleanup_in)
    !> The error variable to be set.
    TYPE(Error_t), INTENT(inout)  :: this
    !> Some string describing the details of the error.
    CHARACTER(len=*), INTENT(in)  :: error_description
    !> If true, the cleanup error handler is called.
    LOGICAL, INTENT(in), OPTIONAL :: immediate_cleanup_in
    !! Local Data
    LOGICAL :: immediate_cleanup

    !! Process Optional Arguments
    immediate_cleanup = .FALSE.
    IF (PRESENT(immediate_cleanup_in)) immediate_cleanup = immediate_cleanup_in

    !! Set Flags and Variables
    this%error_description = error_description
    this%error_set = .TRUE.

    IF (immediate_cleanup) THEN
       CALL Cleanup(this)
    END IF
  END SUBROUTINE SetGenericError
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE SetCustomError(this, error_code, error_description, &
       & immediate_cleanup_in)
    TYPE(Error_t), INTENT(INOUT)  :: this
    INTEGER, INTENT(IN)           :: error_code
    CHARACTER(len=*), INTENT(IN)  :: error_description
    LOGICAL, INTENT(IN), OPTIONAL :: immediate_cleanup_in
    !! Local Data
    LOGICAL :: immediate_cleanup

    immediate_cleanup = .FALSE.
    IF (PRESENT(immediate_cleanup_in)) immediate_cleanup = immediate_cleanup_in

    this%error_description = error_description
    this%error_set = .TRUE.
    this%custom_error = error_code
    this%custom_error_set = .TRUE.

    IF (immediate_cleanup) THEN
       CALL Cleanup(this)
    END IF
  END SUBROUTINE SetCustomError
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Routine to call to check if an MPI error has occurred.
  FUNCTION CheckMPIError(this, error_description, mpi_error, &
       & immediate_cleanup_in) RESULT(error_occurred)
    !> The error variable to be set.
    TYPE(Error_t), INTENT(inout)  :: this
    !> Some string describing the details of the error.
    CHARACTER(len=*), INTENT(in)  :: error_description
    !> The error variable produced by mpi.
    INTEGER, INTENT(in)           :: mpi_error
    !> If true, the cleanup error handler is called.
    LOGICAL, INTENT(in), OPTIONAL :: immediate_cleanup_in
    !> True if an error has occurred, false otherwise.
    LOGICAL :: error_occurred
    !! Local Data
    LOGICAL :: immediate_cleanup

    !! Process Optional Arguments
    immediate_cleanup = .FALSE.
    IF (PRESENT(immediate_cleanup_in)) immediate_cleanup = immediate_cleanup_in

    !! Check Error
    IF (.NOT. mpi_error .EQ. MPI_SUCCESS) THEN
       this%mpi_error_set = .TRUE.
       this%mpi_error = mpi_error
       CALL SetGenericError(this, error_description, immediate_cleanup)
    END IF
    error_occurred = ErrorOccurred(this)
  END FUNCTION CheckMPIError
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Routine to call if an alloc error has occurred.
  FUNCTION CheckAllocError(this, error_description, alloc_error, &
       & immediate_cleanup_in) RESULT(error_occurred)
    !> This the error variable to be set.
    TYPE(Error_t), INTENT(inout)  :: this
    !> Some string describing the details of the error.
    CHARACTER(len=*), INTENT(in)  :: error_description
    !> The error variable produced by alloc.
    INTEGER, INTENT(in)           :: alloc_error
    !> If true, the cleanup error handler is called.
    LOGICAL, INTENT(in), OPTIONAL :: immediate_cleanup_in
    !> True if an error has occurred, false otherwise.
    LOGICAL :: error_occurred
    !! Local Data
    LOGICAL :: immediate_cleanup

    !! Process Optional Arguments
    IF (PRESENT(immediate_cleanup_in)) immediate_cleanup = immediate_cleanup_in

    !! Check Error
    IF (.NOT. alloc_error .EQ. 0) THEN
       this%alloc_error_set = .TRUE.
       this%alloc_error = alloc_error
       CALL SetGenericError(this, error_description, immediate_cleanup)
    END IF
    error_occurred = ErrorOccurred(this)
  END FUNCTION CheckAllocError
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Check if an error has occurred or not.
  FUNCTION ErrorOccurred(this) RESULT(occurred)
    !> The error variable to check.
    TYPE(Error_t), INTENT(in) :: this
    !> True if an error has occurred, false otherwise.
    LOGICAL :: occurred

    occurred = this%error_set
  END FUNCTION ErrorOccurred
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out that an error has occurred.
  SUBROUTINE PrintError(this)
    !> The error to print out.
    TYPE(Error_t), INTENT(in) :: this
    !! Local Data
    CHARACTER(len=80) :: error_string
    INTEGER :: error_string_len
    INTEGER :: error_string_error

    !! Print Out Information About The Error
    IF (ErrorOccurred(this)) THEN
       WRITE(*,'(A)') "# An error has occurred."
       IF (this%alloc_error_set) THEN
          WRITE(*,'(A)') "# Of type: alloc error."
          WRITE(*,'(I3)') this%alloc_error
       ELSE IF (this%mpi_error_set) THEN
          WRITE(*,'(A)') "# Of type: mpi error."
          CALL MPI_Error_String(this%mpi_error, error_string, &
               & error_string_len, error_string_error)
          WRITE(*,'(A)') TRIM(error_string)
       ELSE
          WRITE(*,'(A)') "# Of type: generic error."
       END IF
       WRITE(*,'(A)') "# Details:"
       WRITE(*,'(A)',ADVANCE='no') "#"
       WRITE(*,'(A)') TRIM(this%error_description)
    ELSE
       WRITE(*,'(A)') "# No Error Occured, but PrintError was Called"
    END IF
  END SUBROUTINE PrintError
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> As a last case resort, this will print an error message and quit.
  SUBROUTINE Cleanup(this)
    !> The error which has caused the need to cleanup the program.
    TYPE(Error_t), INTENT(in) :: this
    !! Local Data
    INTEGER :: abort_error

    CALL PrintError(this)
    IF (this%mpi_error_set) THEN
       CALL MPI_Abort(MPI_COMM_WORLD, this%mpi_error, abort_error)
    ELSE
       CALL MPI_Abort(MPI_COMM_WORLD, MPI_ERR_UNKNOWN, abort_error)
    END IF
  END SUBROUTINE Cleanup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE ErrorModule
