!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for writing data to the log file.
MODULE LoggingModule
  USE DataTypesModule, ONLY : NTReal
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER :: CurrentLevel = 0
  LOGICAL :: IsActive = .FALSE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ActivateLogger
  PUBLIC :: DeactivateLogger
  PUBLIC :: EnterSubLog
  PUBLIC :: WriteHeader
  PUBLIC :: WriteListElement
  PUBLIC :: WriteElement
  PUBLIC :: WriteCitation
  PUBLIC :: ExitSubLog
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
  SUBROUTINE ActivateLogger
    IsActive = .TRUE.
  END SUBROUTINE ActivateLogger
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Deactivate the logger.
  SUBROUTINE DeactivateLogger
    IsActive = .FALSE.
  END SUBROUTINE DeactivateLogger
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
  !> Write out a header to the log.
  SUBROUTINE WriteHeader(header_value)
    !> The text of the header.
    CHARACTER(LEN=*), INTENT(IN) :: header_value

    IF (IsActive) THEN
       CALL WriteIndent
       WRITE(*,'(A)',ADVANCE='no') header_value
       WRITE(*,'(A1)') ":"
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

       WRITE(*,'(A)',ADVANCE='no') key
       IF (VALUE) THEN
          WRITE(*,'(A)',ADVANCE='no') ": True"
       ELSE
          WRITE(*,'(A)',ADVANCE='no') ": False"
       END IF

       WRITE(*,*)
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

       WRITE(*,'(A)',ADVANCE='no') key
       WRITE(*,'(A)',ADVANCE='no') ": "
       WRITE(*,'(ES22.14)',ADVANCE='no') VALUE

       WRITE(*,*)
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

       WRITE(*,'(A)',ADVANCE='no') key
       WRITE(*,'(A)',ADVANCE='no') ": "
       WRITE(*,'(I10)',ADVANCE='no') VALUE

       WRITE(*,*)
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

       WRITE(*,'(A)',ADVANCE='no') key
       WRITE(*,'(A)',ADVANCE='no') ": "
       WRITE(*,'(A)',ADVANCE='no') VALUE

       WRITE(*,*)
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

       WRITE(*,'(A)',ADVANCE='no') "- "
       WRITE(*,'(A)',ADVANCE='no') key
       IF (VALUE) THEN
          WRITE(*,'(A)',ADVANCE='no') ": True"
       ELSE
          WRITE(*,'(A)',ADVANCE='no') ": False"
       END IF

       WRITE(*,*)
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

       WRITE(*,'(A)',ADVANCE='no') "- "
       WRITE(*,'(A)',ADVANCE='no') key
       WRITE(*,'(A)',ADVANCE='no') ": "
       WRITE(*,'(ES22.14)',ADVANCE='no') VALUE

       WRITE(*,*)
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

       WRITE(*,'(A)',ADVANCE='no') "- "
       WRITE(*,'(A)',ADVANCE='no') key
       WRITE(*,'(A)',ADVANCE='no') ": "
       WRITE(*,'(I10)',ADVANCE='no') VALUE

       WRITE(*,*)
    END IF
  END SUBROUTINE WriteListElement_int
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write out a list element.
  SUBROUTINE WriteListElement_string(key, VALUE)
    !> Some text to write.
    CHARACTER(LEN=*), INTENT(IN) :: key
    !> A text value to write.
    CHARACTER(LEN=*), INTENT(IN) :: VALUE

    IF (IsActive) THEN
       CALL WriteIndent

       WRITE(*,'(A)',ADVANCE='no') "- "
       WRITE(*,'(A)',ADVANCE='no') key
       WRITE(*,'(A)',ADVANCE='no') ": "
       WRITE(*,'(A)',ADVANCE='no') VALUE

       WRITE(*,*)
    END IF
  END SUBROUTINE WriteListElement_string
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write out a citation element.
  SUBROUTINE WriteCitation(citation_list)
    !> A list of citations, separated by a space.
    CHARACTER(LEN=*), INTENT(IN) :: citation_list
    INTEGER :: pos1, pos2

    IF (IsActive) THEN
       CALL WriteIndent
       WRITE(*,'(A)') "Citations:"
       CALL EnterSubLog

       pos1 = 1
       pos2 = INDEX(citation_list(pos1:), ' ')
       DO WHILE(pos2 .NE. 0)
          CALL WriteIndent
          WRITE(*,'(A)') citation_list(pos1:pos1+pos2-1)
          pos1 = pos1 + pos2
          pos2 = INDEX(citation_list(pos1:), ' ')
       END DO
       CALL WriteIndent
       WRITE(*,'(A)') citation_list(pos1:)

       CALL ExitSubLog
    END IF
  END SUBROUTINE WriteCitation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Writes out the indentation needed for this level
  SUBROUTINE WriteIndent
    INTEGER :: counter

    DO counter=1,CurrentLevel*2
       WRITE(*,'(A1)',ADVANCE='NO') " "
    END DO
  END SUBROUTINE WriteIndent
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE LoggingModule
