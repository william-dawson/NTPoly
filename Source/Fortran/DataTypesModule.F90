!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module to store specifications for basic data types.
MODULE DataTypesModule
  USE NTMPIModule
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_DOUBLE, C_DOUBLE_COMPLEX, C_LONG, &
    & C_FLOAT
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The precision of floating point numbers we will use in this program.
  INTEGER, PARAMETER, PUBLIC :: NTREAL = C_DOUBLE
  !> A low precision type for mixed-precision application.
  INTEGER, PARAMETER, PUBLIC :: NTLOWP = C_FLOAT
  !> MPI floating point datatype with the precision we will use in this program.
  INTEGER, PARAMETER, PUBLIC :: MPINTREAL = MPI_DOUBLE_PRECISION
  !~> MPI low precision floating point for mixed-precision application.
  INTEGER, PARAMETER, PUBLIC :: MPILOWP = MPI_REAL
  !> The complex numbers we will use in this program.
  INTEGER, PARAMETER, PUBLIC :: NTCOMPLEX = C_DOUBLE_COMPLEX
  !> MPI complex datatype with the precision we will use in this program.
  INTEGER, PARAMETER, PUBLIC :: MPINTCOMPLEX = MPI_DOUBLE_COMPLEX
  !> A long integer type for when normal ints will not do
  INTEGER, PARAMETER, PUBLIC :: NTLONG = C_LONG
  !> MPI Integer type we will use in this program.
  INTEGER, PARAMETER, PUBLIC :: MPINTINTEGER = MPI_INTEGER
  !> MPI Integer type we will use in this program.
  INTEGER, PARAMETER, PUBLIC :: MPINTLONG = MPI_INTEGER8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE DataTypesModule
