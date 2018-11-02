!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Simplfiying Per Element Operations on Matrices.
MODULE MatrixMapsModule
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, &
       & GetMatrixTripletList, FillMatrixFromTripletList
  USE TripletListModule, ONLY : TripletList_r, TripletList_c, &
       & DestructTripletList
  USE TripletModule, ONLY : Triplet_r, Triplet_c
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: MapMatrix_psr
  PUBLIC :: MapMatrix_psc
  PUBLIC :: MapTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE MapTripletList
     MODULE PROCEDURE :: MapTripletList_r
     MODULE PROCEDURE :: MapTripletList_c
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Given a distributed matrix, apply this procedure to each element (real).
  SUBROUTINE MapMatrix_psr(inmat, outmat, proc)
    !> The matrix to apply the procedure to.
    TYPE(Matrix_ps), INTENT(IN) :: inmat
    !> The matrix where each element has had proc called on it.
    TYPE(Matrix_ps), INTENT(INOUT) :: outmat
    INTERFACE
       !> The procedure to apply to each element.
       SUBROUTINE proc(row, column, val)
         USE DataTypesModule, ONLY : NTREAL
         !> The row value of an element.
         INTEGER, INTENT(INOUT) :: row
         !> The column value of an element.
         INTEGER, INTENT(INOUT) :: column
         !> The actual value of an element.
         REAL(KIND=NTREAL), INTENT(INOUT) :: val
       END SUBROUTINE proc
    END INTERFACE
    !! Local Variables
    TYPE(TripletList_r) :: inlist, outlist

    INCLUDE "map_includes/MapMatrix.f90"
  END SUBROUTINE MapMatrix_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Given a distributed matrix, apply this procedure to each element (complex).
  SUBROUTINE MapMatrix_psc(inmat, outmat, proc)
    !> The matrix to apply the procedure to.
    TYPE(Matrix_ps), INTENT(IN) :: inmat
    !> The matrix where each element has had proc called on it.
    TYPE(Matrix_ps), INTENT(INOUT) :: outmat
    INTERFACE
       !> The procedure to apply to each element.
       SUBROUTINE proc(row, column, val)
         USE DataTypesModule, ONLY : NTCOMPLEX
         !> The row value of an element.
         INTEGER, INTENT(INOUT) :: row
         !> The column value of an element.
         INTEGER, INTENT(INOUT) :: column
         !> The actual value of an element.
         COMPLEX(KIND=NTCOMPLEX), INTENT(INOUT) :: val
       END SUBROUTINE proc
    END INTERFACE
    !! Local Variables
    TYPE(TripletList_c) :: inlist, outlist

    INCLUDE "map_includes/MapMatrix.f90"
  END SUBROUTINE MapMatrix_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Given a triplet list, apply this procedure to each element.
  SUBROUTINE MapTripletList_r(inlist, outlist, proc)
    !> The matrix to apply the procedure to.
    TYPE(TripletList_r), INTENT(IN) :: inlist
    !> The matrix where each element has had proc called on it.
    TYPE(TripletList_r), INTENT(INOUT) :: outlist
    INTERFACE
       !> The procedure to apply to each element.
       SUBROUTINE proc(row, column, val)
         USE DataTypesModule, ONLY : NTREAL
         !> The row value of an element.
         INTEGER, INTENT(INOUT) :: row
         !> The column value of an element.
         INTEGER, INTENT(INOUT) :: column
         !> The actual value of an element.
         REAL(KIND=NTREAL), INTENT(INOUT) :: val
       END SUBROUTINE proc
    END INTERFACE
    !! Local Variables
    TYPE(Triplet_r) :: temp

    INCLUDE "map_includes/MapTripletList.f90"

  END SUBROUTINE MapTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Given a triplet list, apply this procedure to each element.
  SUBROUTINE MapTripletList_c(inlist, outlist, proc)
    !> The matrix to apply the procedure to.
    TYPE(TripletList_c), INTENT(IN) :: inlist
    !> The matrix where each element has had proc called on it.
    TYPE(TripletList_c), INTENT(INOUT) :: outlist
    INTERFACE
       !> The procedure to apply to each element.
       SUBROUTINE proc(row, column, val)
         USE DataTypesModule, ONLY : NTCOMPLEX
         !> The row value of an element.
         INTEGER, INTENT(INOUT) :: row
         !> The column value of an element.
         INTEGER, INTENT(INOUT) :: column
         !> The actual value of an element.
         COMPLEX(KIND=NTCOMPLEX), INTENT(INOUT) :: val
       END SUBROUTINE proc
    END INTERFACE
    !! Local Variables
    TYPE(Triplet_c) :: temp

    INCLUDE "map_includes/MapTripletList.f90"
  END SUBROUTINE MapTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE MatrixMapsModule
