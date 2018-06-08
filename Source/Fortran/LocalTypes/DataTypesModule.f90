!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module to store specifications for basic data types.
MODULE DataTypesModule
  USE ISO_C_BINDING, ONLY : C_DOUBLE, C_DOUBLE_COMPLEX
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: IsSmallR
  PUBLIC :: IsSmallC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The precision of floating point numbers we will use in this program.
  INTEGER, PARAMETER, PUBLIC :: NTREAL = C_DOUBLE
  !> MPI floating point datatype with the precision we will use in this program.
  INTEGER, PUBLIC :: MPINTREAL = MPI_DOUBLE_PRECISION
  !> The complex numbers we will use in this program.
  INTEGER, PARAMETER, PUBLIC :: NTCOMPLEX = C_DOUBLE_COMPLEX
  !> MPI complex datatype with the precision we will use in this program.
  INTEGER, PARAMETER, PUBLIC :: MPINTCOMPLEX = MPI_DOUBLE_COMPLEX
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Determines if a real variable is small
  PURE FUNCTION IsSmallR(val, thresh) RESULT(issmall)
    !! Parameters
    REAL(NTREAL), INTENT(IN) :: val
    REAL(NTREAL), INTENT(IN) :: thresh
    LOGICAL :: issmall

    issmall = ABS(val) < thresh
  END FUNCTION IsSmallR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Determines if a complex variable is small
  PURE FUNCTION IsSmallC(val, thresh) RESULT(issmall)
    !! Parameters
    COMPLEX(NTCOMPLEX), INTENT(IN) :: val
    REAL(NTREAL), INTENT(IN) :: thresh
    LOGICAL :: issmall
    !! Local

    issmall = ABS(val) < thresh
  END FUNCTION IsSmallC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE DataTypesModule
