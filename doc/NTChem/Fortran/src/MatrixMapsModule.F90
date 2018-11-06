!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Simplfiying Per Element Operations on Matrices.
MODULE MatrixMapsModule
  USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, &
       & GetMatrixTripletList, FillMatrixFromTripletList, &
       & DestructMatrix
  USE TripletListModule, ONLY : TripletList_r, TripletList_c, &
       & DestructTripletList, AppendToTripletList, GetTripletAt, &
       & ConstructTripletList
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
  SUBROUTINE MapMatrix_psr(inmat, outmat, proc, supp_in)
    !> The matrix to apply the procedure to.
    TYPE(Matrix_ps), INTENT(IN) :: inmat
    !> The matrix where each element has had proc called on it.
    TYPE(Matrix_ps), INTENT(INOUT) :: outmat
    INTERFACE
       !> The procedure to apply to each element.
       FUNCTION proc(row, column, val, supp_in) RESULT(valid)
         USE DataTypesModule, ONLY : NTREAL
         !> The row value of an element.
         INTEGER, INTENT(INOUT) :: row
         !> The column value of an element.
         INTEGER, INTENT(INOUT) :: column
         !> The actual value of an element.
         REAL(KIND=NTREAL), INTENT(INOUT) :: val
         !> Any supplementary data you need to pass the map can packed here.
         REAL(KIND=NTREAL), DIMENSION(:), INTENT(IN), OPTIONAL :: supp_in
         !> Set this to false to filter an element.
         LOGICAL :: valid
       END FUNCTION proc
    END INTERFACE
    !> Any supplementary data you need to pass the map can packed here.
    REAL(KIND=NTREAL), DIMENSION(:), INTENT(IN), OPTIONAL :: supp_in
    !! Local Variables
    TYPE(TripletList_r) :: inlist, outlist

    INCLUDE "map_includes/MapMatrix.f90"
  END SUBROUTINE MapMatrix_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Given a distributed matrix, apply this procedure to each element (complex).
  SUBROUTINE MapMatrix_psc(inmat, outmat, proc, supp_in)
    !> The matrix to apply the procedure to.
    TYPE(Matrix_ps), INTENT(IN) :: inmat
    !> The matrix where each element has had proc called on it.
    TYPE(Matrix_ps), INTENT(INOUT) :: outmat
    INTERFACE
       !> The procedure to apply to each element.
       FUNCTION proc(row, column, val, supp_in) RESULT(valid)
         USE DataTypesModule, ONLY : NTCOMPLEX
         !> The row value of an element.
         INTEGER, INTENT(INOUT) :: row
         !> The column value of an element.
         INTEGER, INTENT(INOUT) :: column
         !> The actual value of an element.
         COMPLEX(KIND=NTCOMPLEX), INTENT(INOUT) :: val
         !> Any supplementary data you need to pass the map can packed here.
         COMPLEX(KIND=NTCOMPLEX), DIMENSION(:), INTENT(IN), OPTIONAL :: supp_in
         !> Set this to false to filter an element.
         LOGICAL :: valid
       END FUNCTION proc
    END INTERFACE
    !> Any supplementary data you need to pass the map can packed here.
    COMPLEX(KIND=NTCOMPLEX), DIMENSION(:), INTENT(IN), OPTIONAL :: supp_in
    !! Local Variables
    TYPE(TripletList_c) :: inlist, outlist

    INCLUDE "map_includes/MapMatrix.f90"
  END SUBROUTINE MapMatrix_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Given a triplet list, apply this procedure to each element.
  SUBROUTINE MapTripletList_r(inlist, outlist, proc, supp_in, num_slices_in, &
       & my_slice_in)
    !> The matrix to apply the procedure to.
    TYPE(TripletList_r), INTENT(IN) :: inlist
    !> The matrix where each element has had proc called on it.
    TYPE(TripletList_r), INTENT(INOUT) :: outlist
    INTERFACE
       !> The procedure to apply to each element.
       FUNCTION proc(row, column, val, supp_in) RESULT(valid)
         USE DataTypesModule, ONLY : NTREAL
         !> The row value of an element.
         INTEGER, INTENT(INOUT) :: row
         !> The column value of an element.
         INTEGER, INTENT(INOUT) :: column
         !> The actual value of an element.
         REAL(KIND=NTREAL), INTENT(INOUT) :: val
         !> Any supplementary data you need to pass the map can packed here.
         REAL(KIND=NTREAL), DIMENSION(:), INTENT(IN), OPTIONAL :: supp_in
         !> Set this to false to filter an element.
         LOGICAL :: valid
       END FUNCTION proc
    END INTERFACE
    !> Any supplementary data you need to pass the map can packed here.
    REAL(KIND=NTREAL), DIMENSION(:), INTENT(IN), OPTIONAL :: supp_in
    !> How many process slices to do this mapping on (default is 1)
    INTEGER, INTENT(IN), OPTIONAL :: num_slices_in
    !> What process slice this process should compute (default is 0).
    INTEGER, INTENT(IN), OPTIONAL :: my_slice_in
    !! Local Variables
    TYPE(Triplet_r) :: temp

    INCLUDE "map_includes/MapTripletList.f90"

  END SUBROUTINE MapTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Given a triplet list, apply this procedure to each element.
  SUBROUTINE MapTripletList_c(inlist, outlist, proc, supp_in, num_slices_in, &
       & my_slice_in)
    !> The matrix to apply the procedure to.
    TYPE(TripletList_c), INTENT(IN) :: inlist
    !> The matrix where each element has had proc called on it.
    TYPE(TripletList_c), INTENT(INOUT) :: outlist
    INTERFACE
       !> The procedure to apply to each element.
       FUNCTION proc(row, column, val, supp_in) RESULT(valid)
         USE DataTypesModule, ONLY : NTCOMPLEX
         !> The row value of an element.
         INTEGER, INTENT(INOUT) :: row
         !> The column value of an element.
         INTEGER, INTENT(INOUT) :: column
         !> The actual value of an element.
         COMPLEX(KIND=NTCOMPLEX), INTENT(INOUT) :: val
         !> Any supplementary data you need to pass the map can packed here.
         COMPLEX(KIND=NTCOMPLEX), DIMENSION(:), INTENT(IN), OPTIONAL :: supp_in
         !> Set this to false to filter an element.
         LOGICAL :: valid
       END FUNCTION proc
    END INTERFACE
    !> Any supplementary data you need to pass the map can packed here.
    COMPLEX(KIND=NTCOMPLEX), DIMENSION(:), INTENT(IN), OPTIONAL :: supp_in
    !> How many process slices to do this mapping on (default is 1)
    INTEGER, INTENT(IN), OPTIONAL :: num_slices_in
    !> What process slice this process should compute (default is 0).
    INTEGER, INTENT(IN), OPTIONAL :: my_slice_in
    !! Local Variables
    TYPE(Triplet_c) :: temp

    INCLUDE "map_includes/MapTripletList.f90"
  END SUBROUTINE MapTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE MatrixMapsModule
