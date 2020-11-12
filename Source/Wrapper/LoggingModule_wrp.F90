!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for writing data to the log file.
MODULE LoggingModule_wrp
  USE LoggingModule
  USE ISO_C_BINDING, ONLY : c_char, c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ActivateLogger_wrp
  PUBLIC :: ActivateLoggerFile_wrp
  PUBLIC :: DeactivateLogger_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Activate the logger and have it write to a file.
  SUBROUTINE ActivateLoggerFile_wrp(file_name, name_size) &
       & BIND(C,NAME="ActivateLoggerFile_wrp")
    CHARACTER(KIND=C_CHAR), INTENT(IN) :: file_name(name_size)
    INTEGER(KIND=C_INT), INTENT(IN) :: name_size
    !! Local Data
    CHARACTER(LEN=name_size) :: local_string
    INTEGER :: II

    DO II=1,name_size
       local_string(II:II) = file_name(II)
    END DO

    CALL ActivateLogger(local_string)
  END SUBROUTINE ActivateLoggerFile_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Activate the logger.
  SUBROUTINE ActivateLogger_wrp() BIND(C,NAME="ActivateLogger_wrp")
    CALL ActivateLogger()
  END SUBROUTINE ActivateLogger_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Activate the logger.
  SUBROUTINE DeactivateLogger_wrp() BIND(C,NAME="DeactivateLogger_wrp")
    CALL DeactivateLogger()
  END SUBROUTINE DeactivateLogger_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE LoggingModule_wrp
