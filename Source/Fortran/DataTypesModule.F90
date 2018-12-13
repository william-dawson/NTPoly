!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module to store specifications for basic data types.
MODULE DataTypesModule
  USE NTMPIModule
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_DOUBLE, C_DOUBLE_COMPLEX
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The precision of floating point numbers we will use in this program.
  INTEGER, PARAMETER, PUBLIC :: NTREAL = C_DOUBLE
  !> MPI floating point datatype with the precision we will use in this program.
  INTEGER, PUBLIC :: MPINTREAL = MPI_DOUBLE_PRECISION
  !> The complex numbers we will use in this program.
  INTEGER, PARAMETER, PUBLIC :: NTCOMPLEX = C_DOUBLE_COMPLEX
  !> MPI complex datatype with the precision we will use in this program.
  INTEGER, PARAMETER, PUBLIC :: MPINTCOMPLEX = MPI_DOUBLE_COMPLEX
  !> MPI Integer type we will use in this program.
  INTEGER, PARAMETER, PUBLIC :: MPINTINTEGER = MPI_INT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE DataTypesModule
