!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This module contains helpers for processing matrix market files.
MODULE MatrixMarketModule
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ENUM, BIND(c)
    !> Sparse coordinate file.
    ENUMERATOR :: MM_COORDINATE=1
    !> Dense array file.
    ENUMERATOR :: MM_ARRAY=2
    !> Real data being read in.
    ENUMERATOR :: MM_REAL=1
    !> Integer data being read in.
    ENUMERATOR :: MM_INTEGER=2
    !>Complex numbers being read in.
    ENUMERATOR :: MM_COMPLEX=3
    !> Just a pattern of non zeros.
    ENUMERATOR :: MM_PATTERN=4
    !> File lacks symmetry.
    ENUMERATOR :: MM_GENERAL=1
    !> File is symmetric
    ENUMERATOR :: MM_SYMMETRIC=2
    !> File is skew symmetric.
    ENUMERATOR :: MM_SKEW_SYMMETRIC=3
    !> File is hermitian.
    ENUMERATOR :: MM_HERMITIAN=4
  END ENUM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ParseMMHeader
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Parse a matrix market header.
  !! @param[in]  line
  !! @param[out] sparsity_type sets if coordinate or array type.
  !! @param[out] data_type sets if real, integer, complex, pattern.
  !! @param[out] pattern_type sets if general, symmetric, skew_symmetric,
  !! hermitian.
  !! @result returns true if no errors.
  FUNCTION ParseMMHeader(line,sparsity_type,data_type,pattern_type) &
       & RESULT(no_error)
    !! Parameters
    CHARACTER(len=*), INTENT(IN) :: line
    INTEGER, INTENT(OUT) :: sparsity_type
    INTEGER, INTENT(OUT) :: data_type
    INTEGER, INTENT(OUT) :: pattern_type
    LOGICAL :: no_error
    !! Local Data
    INTEGER :: pos1, pos2

    no_error = .TRUE.

    !! This part is just "MatrixMarket".
    pos1 = 1
    pos2 = INDEX(line(pos1:), ' ')

    !! This part is just "matrix".
    pos1 = pos2+pos1
    pos2 = INDEX(line(pos1:), ' ')

    !! This part is coordinate or array.
    pos1 = pos2+pos1
    pos2 = INDEX(line(pos1:), ' ')
    SELECT CASE(TRIM(line(pos1:pos1+pos2-1)))
    CASE('coordinate')
       sparsity_type = MM_COORDINATE
    CASE('array')
       sparsity_type = MM_ARRAY
    CASE DEFAULT
       no_error = .FALSE.
    END SELECT

    !! This part is real, integer, complex, pattern.
    pos1 = pos2+pos1
    pos2 = INDEX(line(pos1:), ' ')
    SELECT CASE(TRIM(line(pos1:pos1+pos2-1)))
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
END MODULE MatrixMarketModule
