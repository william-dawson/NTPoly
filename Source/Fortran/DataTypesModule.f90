!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module to store specifications for basic data types.
MODULE DataTypesModule
  USE ErrorModule, ONLY : CheckMPIError, Error_t
  USE ISO_C_BINDING
  USE ISO_FORTRAN_ENV
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER, PARAMETER :: BITSPERDOUBLE = 8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The precision of floating point numbers we will use in this program.
  INTEGER, PARAMETER, PUBLIC :: NTREAL = c_double
  !> MPI floating point datatype with the precision we will use in this program.
  INTEGER, PUBLIC :: MPINTREAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: MPITypeInfoInit
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
