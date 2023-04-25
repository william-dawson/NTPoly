!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This module contains helpers for processing matrix market files.
MODULE MatrixMarketModule
  USE DataTypesModule, ONLY : NTREAL, NTLONG
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ENUM, BIND(c)
     !> Sparse coordinate file.
     ENUMERATOR :: MM_COORDINATE = 1
     !> Dense array file.
     ENUMERATOR :: MM_ARRAY = 2
     !> Real data being read in.
     ENUMERATOR :: MM_REAL = 1
     !> Integer data being read in.
     ENUMERATOR :: MM_INTEGER = 2
     !>Complex numbers being read in.
     ENUMERATOR :: MM_COMPLEX = 3
     !> Just a pattern of non zeros.
     ENUMERATOR :: MM_PATTERN = 4
     !> File lacks symmetry.
     ENUMERATOR :: MM_GENERAL = 1
     !> File is symmetric
     ENUMERATOR :: MM_SYMMETRIC = 2
     !> File is skew symmetric.
     ENUMERATOR :: MM_SKEW_SYMMETRIC = 3
     !> File is hermitian. 
     ENUMERATOR :: MM_HERMITIAN = 4
  END ENUM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The longest line size possible according to the spec.
  INTEGER, PARAMETER :: MAX_LINE_LENGTH = 1024
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ParseMMHeader
  PUBLIC :: WriteMMSize
  PUBLIC :: WriteMMLine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE WriteMMLine
     MODULE PROCEDURE WriteMMLine_ii
     MODULE PROCEDURE WriteMMLine_iif
     MODULE PROCEDURE WriteMMLine_iiff
     MODULE PROCEDURE WriteMMLine_f
     MODULE PROCEDURE WriteMMLine_ff
  END INTERFACE WriteMMLine
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Parse a matrix market header.
  FUNCTION ParseMMHeader(line,sparsity_type,data_type,pattern_type) &
       & RESULT(no_error)
    !> String to parse.
    CHARACTER(LEN = *), INTENT(IN) :: line
    !> If coordinate or array type.
    INTEGER, INTENT(OUT) :: sparsity_type
    !> If real, integer, complex, pattern.
    INTEGER, INTENT(OUT) :: data_type
    !> If general, symmetric, skew_symmetric, hermitian.
    INTEGER, INTENT(OUT) :: pattern_type
    !> True if no errors.
    LOGICAL :: no_error
    !! Local Data
    INTEGER :: pos1, pos2

    no_error = .TRUE.

    !! This part is just "MatrixMarket".
    pos1 = 1
    pos2 = INDEX(line(pos1:), ' ')

    !! This part is just "matrix".
    pos1 = pos2 + pos1
    pos2 = INDEX(line(pos1:), ' ')

    !! This part is coordinate or array.
    pos1 = pos2 + pos1
    pos2 = INDEX(line(pos1:), ' ')
    SELECT CASE(TRIM(line(pos1:pos1 + pos2 - 1)))
    CASE('coordinate')
       sparsity_type = MM_COORDINATE
    CASE('array')
       sparsity_type = MM_ARRAY
    CASE DEFAULT
       no_error = .FALSE.
    END SELECT

    !! This part is real, integer, complex, pattern.
    pos1 = pos2 + pos1
    pos2 = INDEX(line(pos1:), ' ')
    SELECT CASE(TRIM(line(pos1: pos1 + pos2 - 1)))
    CASE('real')
       data_type = MM_REAL
    CASE('array')
       data_type = MM_INTEGER
    CASE('complex')
       data_type = MM_COMPLEX
    CASE('pattern')
       data_type = MM_PATTERN
    CASE DEFAULT
       no_error = .FALSE.
    END SELECT

    !! This part is general, symmetric, skew-symmetric, hermitian.
    pos1 = pos2+pos1
    SELECT CASE(TRIM(line(pos1:)))
    CASE('general')
       pattern_type = MM_GENERAL
    CASE('symmetric')
       pattern_type = MM_SYMMETRIC
    CASE('skew-symmetric')
       pattern_type = MM_SKEW_SYMMETRIC
    CASE('hermitian')
       pattern_type = MM_HERMITIAN
    CASE DEFAULT
       no_error = .FALSE.
    END SELECT

  END FUNCTION ParseMMHeader
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write the line describing the size of the matrix
  PURE SUBROUTINE WriteMMSize(outstring, rows, columns, values_in)
    !> The final string is written to this variable.
    CHARACTER(LEN = MAX_LINE_LENGTH), INTENT(INOUT) :: outstring
    !> The number of rows of the matrix
    INTEGER, INTENT(IN) :: rows
    !> The number of columns of the matrix
    INTEGER, INTENT(IN) :: columns
    !> The total number of non zero values in the matrix (for sparse format).
    INTEGER(KIND = NTLONG), INTENT(IN), OPTIONAL :: values_in
    !! Local variables
    CHARACTER(LEN = MAX_LINE_LENGTH) :: temp1, temp2, temp3

    !! Write everything to strings.
    WRITE(temp1, *) rows
    WRITE(temp2, *) columns
    IF (PRESENT(values_in)) THEN
       WRITE(temp3, *) values_in
    ELSE
       WRITE(temp3, *) ""
    END IF

    !! Combine
    WRITE(outstring, '(3A)') ADJUSTL(TRIM(temp1)), ADJUSTL(TRIM(temp2)), &
         & ADJUSTL(TRIM(temp3))

  END SUBROUTINE WriteMMSize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write a single line that would correspond to a matrix market entry.
  PURE SUBROUTINE WriteMMLine_ii(outstring, row, column, add_newline_in)
    !> The final string is written to this variable.
    CHARACTER(LEN = MAX_LINE_LENGTH), INTENT(INOUT) :: outstring
    !> The first coordinate value
    INTEGER, INTENT(IN) :: row
    !> The second coordinate value
    INTEGER, INTENT(IN) :: column
    !> Whether to append a new line to the output (default = .false.)
    LOGICAL, INTENT(IN), OPTIONAL :: add_newline_in
    !! Local variables
    CHARACTER(LEN = MAX_LINE_LENGTH) :: temp1, temp2
    LOGICAL :: add_newline

    !! Process Optional Arguments
    IF (PRESENT(add_newline_in)) THEN
       add_newline = add_newline_in
    ELSE
       add_newline = .FALSE.
    END IF

    !! Write everything to strings.
    WRITE(temp1, *) row
    WRITE(temp2, *) column

    !! Combine
    IF (add_newline) THEN
       WRITE(outstring, '(3A)') ADJUSTL(TRIM(temp1)), &
            & ADJUSTL(TRIM(temp2)) // NEW_LINE('A')
    ELSE
       WRITE(outstring, '(2A)') ADJUSTL(TRIM(temp1)), ADJUSTL(TRIM(temp2))
    END IF
  END SUBROUTINE WriteMMLine_ii
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write a single line that would correspond to a matrix market entry.
  PURE SUBROUTINE WriteMMLine_iif(outstring, row, column, val, add_newline_in)
    !> The final string is written to this variable.
    CHARACTER(LEN = MAX_LINE_LENGTH), INTENT(INOUT) :: outstring
    !> The first coordinate value
    INTEGER, INTENT(IN) :: row
    !> The second coordinate value
    INTEGER, INTENT(IN) :: column
    !> The value at that coordinate
    REAL(NTREAL), INTENT(IN) :: val
    !> Whether to append a new line to the output (default = .false.)
    LOGICAL, INTENT(IN), OPTIONAL :: add_newline_in
    !! Local variables
    CHARACTER(LEN = MAX_LINE_LENGTH) :: temp1, temp2, temp3
    LOGICAL :: add_newline

    !! Process Optional Arguments
    IF (PRESENT(add_newline_in)) THEN
       add_newline = add_newline_in
    ELSE
       add_newline = .FALSE.
    END IF

    !! Write everything to strings.
    WRITE(temp1, *) row
    WRITE(temp2, *) column
    WRITE(temp3, *) val

    !! Combine
    IF (add_newline) THEN
       WRITE(outstring, '(4A)') ADJUSTL(TRIM(temp1)), &
            & ADJUSTL(TRIM(temp2)), ADJUSTL(TRIM(temp3)) // NEW_LINE('A')
    ELSE
       WRITE(outstring, '(3A)') ADJUSTL(TRIM(temp1)), ADJUSTL(TRIM(temp2)), &
            & ADJUSTL(TRIM(temp3))
    END IF
  END SUBROUTINE WriteMMLine_iif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write a single line that would correspond to a matrix market entry.
  PURE SUBROUTINE WriteMMLine_iiff(outstring, row, column, val1, val2, &
       & add_newline_in)
    !> The final string is written to this variable.
    CHARACTER(LEN = MAX_LINE_LENGTH), INTENT(INOUT) :: outstring
    !> The first coordinate value
    INTEGER, INTENT(IN) :: row
    !> The second coordinate value
    INTEGER, INTENT(IN) :: column
    !> The value at that coordinate
    REAL(NTREAL), INTENT(IN) :: val1
    !> The second value at the coordinate
    REAL(NTREAL), INTENT(IN) :: val2
    !> Whether to append a new line to the output (default = .false.)
    LOGICAL, INTENT(IN), OPTIONAL :: add_newline_in
    !! Local variables
    CHARACTER(LEN = MAX_LINE_LENGTH) :: temp1, temp2, temp3, temp4
    LOGICAL :: add_newline

    !! Process Optional Arguments
    IF (PRESENT(add_newline_in)) THEN
       add_newline = add_newline_in
    ELSE
       add_newline = .FALSE.
    END IF

    !! Write everything to strings.
    WRITE(temp1, *) row
    WRITE(temp2, *) column
    WRITE(temp3, *) val1
    WRITE(temp4, *) val2

    !! Combine
    IF (add_newline) THEN
       WRITE(outstring, '(5A)') ADJUSTL(TRIM(temp1)), &
            & ADJUSTL(TRIM(temp2)), ADJUSTL(TRIM(temp3)), &
            & ADJUSTL(TRIM(temp4)) // NEW_LINE('A')
    ELSE
       WRITE(outstring, '(4A)') ADJUSTL(TRIM(temp1)), &
            & ADJUSTL(TRIM(temp2)), ADJUSTL(TRIM(temp3)), ADJUSTL(TRIM(temp4))
    END IF
  END SUBROUTINE WriteMMLine_iiff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write a single line that would correspond to a matrix market entry.
  PURE SUBROUTINE WriteMMLine_f(outstring, val, add_newline_in)
    !> The final string is written to this variable.
    CHARACTER(LEN = MAX_LINE_LENGTH), INTENT(INOUT) :: outstring
    !> The value at that coordinate
    REAL(NTREAL), INTENT(IN) :: val
    !> Whether to append a new line to the output (default = .false.)
    LOGICAL, INTENT(IN), OPTIONAL :: add_newline_in
    !! Local Variables
    CHARACTER(LEN = MAX_LINE_LENGTH) :: temp1
    LOGICAL :: add_newline

    !! Process Optional Arguments
    IF (PRESENT(add_newline_in)) THEN
       add_newline = add_newline_in
    ELSE
       add_newline = .FALSE.
    END IF

    !! Write everything to strings.
    WRITE(temp1, *) val

    !! Combine
    IF (add_newline) THEN
       WRITE(outstring, '(2A)') ADJUSTL(TRIM(temp1)) // NEW_LINE('A')
    ELSE
       WRITE(outstring, '(A)') ADJUSTL(TRIM(temp1))
    END IF
  END SUBROUTINE WriteMMLine_f
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write a single line that would correspond to a matrix market entry.
  PURE SUBROUTINE WriteMMLine_ff(outstring, val1, val2, add_newline_in)
    !> The final string is written to this variable.
    CHARACTER(LEN = MAX_LINE_LENGTH), INTENT(INOUT) :: outstring
    !> The value at that coordinate
    REAL(NTREAL), INTENT(IN) :: val1
    !> The second value at that coordinate
    REAL(NTREAL), INTENT(IN) :: val2
    !> Whether to append a new line to the output (default = .false.)
    LOGICAL, INTENT(IN), OPTIONAL :: add_newline_in
    !! Local variables
    CHARACTER(LEN = MAX_LINE_LENGTH) :: temp1, temp2
    LOGICAL :: add_newline

    !! Process Optional Arguments
    IF (PRESENT(add_newline_in)) THEN
       add_newline = add_newline_in
    ELSE
       add_newline = .FALSE.
    END IF

    !! Write everything to strings.
    WRITE(temp1, *) val1
    WRITE(temp2, *) val2

    !! Combine
    IF (add_newline) THEN
       WRITE(outstring, '(3A)') ADJUSTL(TRIM(temp1)), ADJUSTL(TRIM(temp2)) &
            & // NEW_LINE('A')
    ELSE
       WRITE(outstring, '(2A)') ADJUSTL(TRIM(temp1)), ADJUSTL(TRIM(temp2))
    END IF
  END SUBROUTINE WriteMMLine_ff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE MatrixMarketModule
