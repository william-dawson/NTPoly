!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for writing data to the log file.
MODULE LoggingModule
  USE DataTypesModule, ONLY : NTReal
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER, SAVE :: CurrentLevel = 0
  LOGICAL, SAVE :: IsActive = .FALSE.
  INTEGER, SAVE :: UNIT = 6
  LOGICAL, SAVE :: file_open = .FALSE.
  INTEGER, SAVE :: initial_offset = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ActivateLogger
  PUBLIC :: DeactivateLogger
  PUBLIC :: IsLoggerActive
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: EnterSubLog
  PUBLIC :: ExitSubLog
  PUBLIC :: WriteElement
  PUBLIC :: WriteHeader
  PUBLIC :: WriteListElement
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SetInitialOffset
  PUBLIC :: SetLoggerLevel
  PUBLIC :: GetLoggerLevel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE WriteListElement
     MODULE PROCEDURE WriteListElement_bool
     MODULE PROCEDURE WriteListElement_float
     MODULE PROCEDURE WriteListElement_int
     MODULE PROCEDURE WriteListElement_string
  END INTERFACE WriteListElement
  INTERFACE WriteElement
     MODULE PROCEDURE WriteElement_bool
     MODULE PROCEDURE WriteElement_float
     MODULE PROCEDURE WriteElement_int
     MODULE PROCEDURE WriteElement_string
  END INTERFACE WriteElement
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Activate the logger.
  SUBROUTINE ActivateLogger(start_document_in, file_name_in, unit_in)
    !> If this is a new document we can write the start document marker.
    LOGICAL, INTENT(IN), OPTIONAL :: start_document_in
    !> An optional file name for writing to.
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: file_name_in
    !> An optional fortran i/o unit override.
    INTEGER, INTENT(IN), OPTIONAL :: unit_in

    IsActive = .TRUE.

    IF (PRESENT(unit_in)) THEN
       UNIT = unit_in
    END IF

    IF (PRESENT(file_name_in)) THEN
       IF (.NOT. PRESENT(unit_in)) THEN
          UNIT = 14
       END IF
       OPEN(unit = UNIT, file = file_name_in)
       file_open = .TRUE.
    END IF

    IF (PRESENT(start_document_in)) THEN
       IF (start_document_in) THEN
          WRITE(UNIT, '(A3)') "---"
          initial_offset = 1
       END IF
    END IF
  END SUBROUTINE ActivateLogger
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Deactivate the logger.
  SUBROUTINE DeactivateLogger
    IsActive = .FALSE.
    IF (file_open) THEN
       CLOSE(UNIT)
    END IF
    UNIT = 6
    CurrentLevel = 0
  END SUBROUTINE DeactivateLogger
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Check if the logger is currently active
  FUNCTION IsLoggerActive() RESULT(active)
    LOGICAL :: active

    active = IsActive
  END FUNCTION IsLoggerActive
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Call this subroutine when you enter into a section with verbose output
  SUBROUTINE EnterSubLog
    CurrentLevel = CurrentLevel + 1
  END SUBROUTINE EnterSubLog
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Call this subroutine when you exit a section with verbose output
  SUBROUTINE ExitSubLog
    CurrentLevel = CurrentLevel - 1
  END SUBROUTINE ExitSubLog
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set a manual initial offset spacing.
  SUBROUTINE SetInitialOffset(offset)
    !> Number of spaces to offset
    INTEGER, INTENT(IN) :: offset

    initial_offset = offset
  END SUBROUTINE SetInitialOffset
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write out a header to the log.
  SUBROUTINE WriteHeader(header_value)
    !> The text of the header.
    CHARACTER(LEN=*), INTENT(IN) :: header_value

    IF (IsActive) THEN
       CALL WriteIndent
       WRITE(UNIT, '(A)', ADVANCE='no') header_value
       WRITE(UNIT, '(A1)') ":"
    END IF
  END SUBROUTINE WriteHeader
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write out a element.
  SUBROUTINE WriteElement_bool(key, VALUE)
    !> Some text to write.
    CHARACTER(LEN=*), INTENT(IN) :: key
    !> An integer value to write.
    LOGICAL, INTENT(IN) :: VALUE

    IF (IsActive) THEN
       CALL WriteIndent

       WRITE(UNIT, '(A)', ADVANCE='no') key
       IF (VALUE) THEN
          WRITE(UNIT, '(A)', ADVANCE='no') ": True"
       ELSE
          WRITE(UNIT, '(A)', ADVANCE='no') ": False"
       END IF

       WRITE(UNIT,*)
    END IF
  END SUBROUTINE WriteElement_bool
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write out a element.
  SUBROUTINE WriteElement_float(key, VALUE)
    !> Some text to write.
    CHARACTER(LEN=*), INTENT(IN) :: key
    !> A float value to write.
    REAL(NTReal), INTENT(IN) :: VALUE

    IF (IsActive) THEN
       CALL WriteIndent

       WRITE(UNIT, '(A)', ADVANCE='no') key
       WRITE(UNIT, '(A)', ADVANCE='no') ": "
       WRITE(UNIT, '(ES22.14)', ADVANCE='no') VALUE

       WRITE(UNIT, *)
    END IF
  END SUBROUTINE WriteElement_float
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write out a element.
  SUBROUTINE WriteElement_int(key, VALUE)
    !> Some text to write.
    CHARACTER(LEN=*), INTENT(IN) :: key
    !> An integer value to write.
    INTEGER, INTENT(IN) :: VALUE

    IF (IsActive) THEN
       CALL WriteIndent

       WRITE(UNIT, '(A)', ADVANCE = 'no') key
       WRITE(UNIT, '(A)', ADVANCE = 'no') ": "
       WRITE(UNIT, '(I20)', ADVANCE = 'no') VALUE

       WRITE(UNIT, *)
    END IF
  END SUBROUTINE WriteElement_int
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write out a element.
  SUBROUTINE WriteElement_string(key, VALUE)
    !> Some text to write.
    CHARACTER(LEN=*), INTENT(IN) :: key
    !> A text value to write.
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: VALUE

    IF (IsActive) THEN
       CALL WriteIndent

       WRITE(UNIT, '(A)', ADVANCE = 'no') key
       WRITE(UNIT, '(A)', ADVANCE = 'no') ": "
       WRITE(UNIT, '(A)', ADVANCE = 'no') VALUE

       WRITE(UNIT,*)
    END IF
  END SUBROUTINE WriteElement_string
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write out a list element.
  SUBROUTINE WriteListElement_bool(key, VALUE)
    !> Some text to write.
    CHARACTER(LEN=*), INTENT(IN) :: key
    !> A bool value to write.
    LOGICAL, INTENT(IN) :: VALUE

    IF (IsActive) THEN
       CALL WriteIndent

       WRITE(UNIT, '(A)', ADVANCE = 'no') "- "
       WRITE(UNIT, '(A)', ADVANCE = 'no') key
       IF (VALUE) THEN
          WRITE(UNIT, '(A)', ADVANCE = 'no') ": True"
       ELSE
          WRITE(UNIT, '(A)', ADVANCE = 'no') ": False"
       END IF

       WRITE(UNIT, *)
    END IF
  END SUBROUTINE WriteListElement_bool
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write out a list element.
  SUBROUTINE WriteListElement_float(key, VALUE)
    !> Some text to write.
    CHARACTER(LEN=*), INTENT(IN) :: key
    !> A float value to write.
    REAL(NTReal), INTENT(IN) :: VALUE

    IF (IsActive) THEN
       CALL WriteIndent

       WRITE(UNIT, '(A)', ADVANCE = 'no') "- "
       WRITE(UNIT, '(A)', ADVANCE = 'no') key
       WRITE(UNIT, '(A)', ADVANCE = 'no') ": "
       WRITE(UNIT, '(ES22.14)', ADVANCE = 'no') VALUE

       WRITE(UNIT,*)
    END IF
  END SUBROUTINE WriteListElement_float
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write out a list element.
  SUBROUTINE WriteListElement_int(key, VALUE)
    !> Some text to write.
    CHARACTER(LEN=*), INTENT(IN) :: key
    !> An integer value to write.
    INTEGER, INTENT(IN) :: VALUE

    IF (IsActive) THEN
       CALL WriteIndent

       WRITE(UNIT, '(A)', ADVANCE = 'no') "- "
       WRITE(UNIT, '(A)', ADVANCE = 'no') key
       WRITE(UNIT, '(A)', ADVANCE = 'no') ": "
       WRITE(UNIT, '(I10)', ADVANCE = 'no') VALUE

       WRITE(UNIT,*)
    END IF
  END SUBROUTINE WriteListElement_int
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write out a list element.
  SUBROUTINE WriteListElement_string(key, VALUE)
    !> Some text to write.
    CHARACTER(LEN=*), INTENT(IN) :: key
    !> A text value to write.
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: VALUE

    IF (IsActive) THEN
       CALL WriteIndent

       WRITE(UNIT, '(A)', ADVANCE = 'no') "- "
       WRITE(UNIT, '(A)', ADVANCE = 'no') key
       IF (PRESENT(VALUE)) THEN
          WRITE(UNIT, '(A)', ADVANCE = 'no') ": "
          WRITE(UNIT, '(A)', ADVANCE = 'no') VALUE
       END IF

       WRITE(UNIT,*)
    END IF
  END SUBROUTINE WriteListElement_string
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Writes out the indentation needed for this level
  SUBROUTINE WriteIndent
    INTEGER :: II

    DO II=1,initial_offset
       WRITE(UNIT, '(A1)', ADVANCE = 'NO') " "
    END DO
    DO II=1,CurrentLevel*2
       WRITE(UNIT, '(A1)', ADVANCE = 'NO') " "
    END DO
  END SUBROUTINE WriteIndent
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the logging level manually
  SUBROUTINE SetLoggerLevel(level)
    INTEGER, INTENT(IN) :: level

    CurrentLevel = level
  END SUBROUTINE SetLoggerLevel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the current logging level
  FUNCTION GetLoggerLevel() RESULT(level)
    INTEGER :: level

    level = CurrentLevel
  END FUNCTION GetLoggerLevel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE LoggingModule