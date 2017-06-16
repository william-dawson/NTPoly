!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module to store specifications for basic data types.
MODULE DataTypesModule
  USE ErrorModule
  USE iso_c_Binding
  USE ISO_FORTRAN_ENV
  USE mpi
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER, PARAMETER :: BITSPERDOUBLE = 8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The precision of floating point numbers we will use in this program.
  INTEGER, PARAMETER, PUBLIC :: NTREAL = c_double
  !integer, parameter, public :: NTREAL = selected_real_kind( &
  !    & BITSPERDOUBLE)
  !> MPI floating point datatype with the precision we will use in this program.
  INTEGER, PUBLIC :: MPINTREAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: MPITypeInfoInit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A subroutine to make a double precision MPI type with the requested
  !! precision. Must be called before MPINTREAL is used.
  SUBROUTINE MPITypeInfoInit()
    !! Local Data
    TYPE(Error_t) :: error
    INTEGER :: mpi_error
    LOGICAL :: error_occurred

    CALL MPI_TYPE_CREATE_F90_REAL(BITSPERDOUBLE,MPI_UNDEFINED, &
         & MPINTREAL, mpi_error)
    MPINTREAL = MPI_DOUBLE_PRECISION
    error_occurred = CheckMPIError(error,"Creating MPI Data Types",mpi_error, &
         & immediate_cleanup_in = .TRUE.)
  END SUBROUTINE MPITypeInfoInit
END MODULE DataTypesModule
