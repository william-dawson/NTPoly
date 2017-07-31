!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for writing data to the log file.
MODULE LoggingModule
  USE DataTypesModule, ONLY : NTReal
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER :: CurrentLevel=0
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
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Activate the logger.
  SUBROUTINE ActivateLogger
    IsActive = .TRUE.
  END SUBROUTINE ActivateLogger
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Deactivate the logger.
  SUBROUTINE DeactivateLogger
    IsActive = .TRUE.
  END SUBROUTINE DeactivateLogger
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Call this function when you enter into a section with verbose output
  SUBROUTINE EnterSubLog
    CurrentLevel = CurrentLevel + 1
  END SUBROUTINE EnterSubLog
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Call this function when you exit a section with verbose output
  SUBROUTINE ExitSubLog
    CurrentLevel = CurrentLevel - 1
  END SUBROUTINE ExitSubLog
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write out a header to the log.
  !! @param[in] header_value the text of the header.
  SUBROUTINE WriteHeader(header_value)
    !! Parameters
    CHARACTER(LEN=*), INTENT(IN) :: header_value

    IF (IsActive) THEN
      CALL WriteIndent
      WRITE(*,'(A)',ADVANCE='no') header_value
      WRITE(*,'(A1)') ":"
    END IF
  END SUBROUTINE WriteHeader
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write out a element.
  !! @param[in] key some text to write.
  !! @param[in] text_value_in a text value to write.
  !! @param[in] int_value_in an integer value to write.
  !! @param[in] float_value_in an float value to write.
  !! @param[in] bool_value_in a bool value to write.
  SUBROUTINE WriteElement(key, text_value_in, int_value_in, float_value_in, &
    & bool_value_in)
    !! Parameters
    CHARACTER(LEN=*), INTENT(IN) :: key
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: text_value_in
    INTEGER, INTENT(IN), OPTIONAL :: int_value_in
    REAL(NTReal), INTENT(IN), OPTIONAL :: float_value_in
    LOGICAL, INTENT(IN), OPTIONAL :: bool_value_in

    IF (IsActive) THEN
      CALL WriteIndent

      WRITE(*,'(A)',ADVANCE='no') key

      IF (PRESENT(text_value_in)) THEN
         WRITE(*,'(A)',ADVANCE='no') ": "
         WRITE(*,'(A)',ADVANCE='no') text_value_in
      END IF
      IF (PRESENT(int_value_in)) THEN
         WRITE(*,'(A)',ADVANCE='no') ": "
         WRITE(*,'(I10)',ADVANCE='no') int_value_in
      END IF
      IF (PRESENT(float_value_in)) THEN
         WRITE(*,'(A)',ADVANCE='no') ": "
         WRITE(*,'(ES10.3)',ADVANCE='no') float_value_in
      END IF
      IF (PRESENT(bool_value_in)) THEN
         IF (bool_value_in) THEN
           WRITE(*,'(A)',ADVANCE='no') ": True"
         ELSE
           WRITE(*,'(A)',ADVANCE='no') ": False"
         END IF
      END IF

      WRITE(*,*)
    END IF
  END SUBROUTINE WriteElement
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write out a list element.
  !! @param[in] key some text to add to the list.
  !! @param[in] text_value_in a text value to add to the list.
  !! @param[in] int_value_in an integer value to add to the list.
  !! @param[in] float_value_in a float value to add to the list.
  !! @param[in] bool_value_in a bool value to add to the list.
  SUBROUTINE WriteListElement(key, text_value_in, int_value_in, float_value_in,&
    bool_value_in)
    !! Parameters
    CHARACTER(LEN=*), INTENT(IN) :: key
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: text_value_in
    INTEGER, INTENT(IN), OPTIONAL :: int_value_in
    REAL(NTReal), INTENT(IN), OPTIONAL :: float_value_in
    LOGICAL, INTENT(IN), OPTIONAL :: bool_value_in

    IF (IsActive) THEN
      CALL WriteIndent

      WRITE(*,'(A)',ADVANCE='no') "- "
      WRITE(*,'(A)',ADVANCE='no') key

      IF (PRESENT(text_value_in)) THEN
         WRITE(*,'(A)',ADVANCE='no') ": "
         WRITE(*,'(A)',ADVANCE='no') text_value_in
      END IF
      IF (PRESENT(int_value_in)) THEN
         WRITE(*,'(A)',ADVANCE='no') ": "
         WRITE(*,'(I10)',ADVANCE='no') int_value_in
      END IF
      IF (PRESENT(float_value_in)) THEN
         WRITE(*,'(A)',ADVANCE='no') ": "
         WRITE(*,'(ES11.4)',ADVANCE='no') float_value_in
      END IF
      IF (PRESENT(bool_value_in)) THEN
         IF (bool_value_in) THEN
           WRITE(*,'(A)',ADVANCE='no') ": True"
         ELSE
           WRITE(*,'(A)',ADVANCE='no') ": False"
         END IF
      END IF

      WRITE(*,*)
    END IF
  END SUBROUTINE WriteListElement
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write out a citation element.
  !! @param[in] citation_list list of citations, separated by a space.
  SUBROUTINE WriteCitation(citation_list)
    !! Parameters
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
END MODULE LoggingModule
