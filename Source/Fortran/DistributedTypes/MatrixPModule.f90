!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Performing Distributed Matrix Operations.
MODULE MatrixPModule
  USE ProcessGridModule, ONLY : ProcessGrid_t
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE, ABSTRACT, PUBLIC :: Matrix_p
     !> Number of matrix rows/columns for full matrix, scaled for process grid.
     INTEGER :: logical_matrix_dimension
     !> Number of matrix rows/columns for the full matrix, unscaled.
     INTEGER :: actual_matrix_dimension
     !! Local Storage
     INTEGER :: start_column !< first column stored locally.
     INTEGER :: end_column !< last column stored locally  is less than this.
     INTEGER :: start_row !< first row stored locally.
     INTEGER :: end_row !< last row stored locally is less than this.
     INTEGER :: local_columns !< number of local columns.
     INTEGER :: local_rows !< number of local rows.
     !! Process Grid
     TYPE(ProcessGrid_t) :: process_grid !< process grid to operate on
  CONTAINS
     !! Construct/Destruct
     PROCEDURE(ConstructEmptyMatrix_p), DEFERRED :: InitEmpty
     PROCEDURE(DestructMatrix_p), DEFERRED :: Destruct
     PROCEDURE(CopyMatrix_p), DEFERRED :: Copy
     !! File I/O
     PROCEDURE(ConstructFromMatrixMarket_p), DEFERRED :: InitMatrixMarket
     PROCEDURE(ConstructFromBinary_p), DEFERRED :: InitBinary
     PROCEDURE(WriteMatrixToMatrixMarket_p), DEFERRED :: WriteMatrixMarket
     PROCEDURE(WriteMatrixToBinary_p), DEFERRED :: WriteBinary
     !! Special Matrices
     PROCEDURE(FillMatrixFromTripletList_p), DEFERRED :: FillFromTripletList
     PROCEDURE(FillMatrixIdentity_p), DEFERRED :: FillIdentity
     !! Basic Accessors
     PROCEDURE(GetTripletList_p), DEFERRED :: GetTripletList
     PROCEDURE(GetBlock_p), DEFERRED :: GetBlock
     !! Basic Accessors
     PROCEDURE :: ActualDim => GetMatrixActualDimension_p
     PROCEDURE :: LogicalDim => GetMatrixLogicalDimension_p
     !! Printing To The Console
     PROCEDURE(PrintMatrix_p), DEFERRED :: Print
     PROCEDURE(PrintMatrixInformation_p), DEFERRED :: PrintInfo
     !! Utilities
     PROCEDURE(TransposeMatrix_p), DEFERRED :: Transpose
     PROCEDURE(CommSplitMatrix_p), DEFERRED :: CommSplit
  END TYPE Matrix_p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ABSTRACT INTERFACE
     SUBROUTINE ConstructEmptyMatrix_p(this, matrix_dim, process_grid_in)
       USE ProcessGridModule, ONLY : ProcessGrid_t
       IMPORT :: Matrix_p
       IMPLICIT NONE
       CLASS(Matrix_p), INTENT(INOUT) :: this
       INTEGER, INTENT(IN)            :: matrix_dim
       TYPE(ProcessGrid_t), INTENT(IN), OPTIONAL :: process_grid_in
     END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     PURE SUBROUTINE DestructMatrix_p(this)
       IMPORT :: Matrix_p
       IMPLICIT NONE
       CLASS(Matrix_p), INTENT(INOUT) :: this
     END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE CopyMatrix_p(this, matA)
       IMPORT :: Matrix_p
       IMPLICIT NONE
       CLASS(Matrix_p), INTENT(INOUT) :: this
       CLASS(Matrix_p), INTENT(IN)    :: matA
     END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE ConstructFromMatrixMarket_p(this, file_name)
       IMPORT :: Matrix_p
       IMPLICIT NONE
       !! Parameters
       CLASS(Matrix_p), INTENT(INOUT) :: this
       CHARACTER(len=*), INTENT(IN)  :: file_name
     END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE ConstructFromBinary_p(this, file_name)
       IMPORT :: Matrix_p
       IMPLICIT NONE
       !! Parameters
       CLASS(Matrix_p), INTENT(INOUT) :: this
       CHARACTER(len=*), INTENT(IN)  :: file_name
     END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE WriteMatrixToMatrixMarket_p(this, file_name)
       IMPORT :: Matrix_p
       IMPLICIT NONE
       !! Parameters
       CLASS(Matrix_p), INTENT(IN)   :: this
       CHARACTER(len=*), INTENT(IN) :: file_name
     END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE WriteMatrixToBinary_p(this, file_name)
       IMPORT :: Matrix_p
       IMPLICIT NONE
       !! Parameters
       CLASS(Matrix_p), INTENT(IN)   :: this
       CHARACTER(len=*), INTENT(IN) :: file_name
     END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE FillMatrixFromTripletList_p(this, triplet_list, &
       & preduplicated_in)
       USE TripletListModule, ONLY : TripletList
       IMPORT :: Matrix_p
       IMPLICIT NONE
       CLASS(Matrix_p), INTENT(INOUT) :: this
       CLASS(TripletList), INTENT(IN) :: triplet_list
       LOGICAL, INTENT(IN), OPTIONAL :: preduplicated_in
    END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE FillMatrixIdentity_p(this)
       IMPORT :: Matrix_p
       IMPLICIT NONE
       CLASS(Matrix_p), INTENT(INOUT) :: this
     END SUBROUTINE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     PURE SUBROUTINE GetTripletList_p(this, triplet_list)
        USE TripletListModule, ONLY : TripletList
        IMPORT :: Matrix_p
        IMPLICIT NONE
        CLASS(Matrix_p), INTENT(IN) :: this
        CLASS(TripletList), INTENT(INOUT) :: triplet_list
     END SUBROUTINE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE GetBlock_p(this, triplet_list, start_row, end_row, &
        & start_column, end_column)
        USE TripletListModule, ONLY : TripletList
        IMPORT :: Matrix_p
        IMPLICIT NONE
        CLASS(Matrix_p), INTENT(IN) :: this
        CLASS(TripletList), INTENT(INOUT) :: triplet_list
        INTEGER, INTENT(IN) :: start_row, end_row
        INTEGER, INTENT(IN) :: start_column, end_column
     END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE PrintMatrix_p(this, file_name_in)
       IMPORT :: Matrix_p
       IMPLICIT NONE
       CLASS(Matrix_p), INTENT(IN) :: this
       CHARACTER(len=*), OPTIONAL, INTENT(IN) :: file_name_in
     END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE PrintMatrixInformation_p(this)
       IMPORT :: Matrix_p
       IMPLICIT NONE
       CLASS(Matrix_p), INTENT(IN) :: this
     END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE TransposeMatrix_p(this, matA)
       IMPORT :: Matrix_p
       IMPLICIT NONE
       CLASS(Matrix_p), INTENT(INOUT) :: this
       CLASS(Matrix_p), INTENT(IN) :: matA
     END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE CommSplitMatrix_p(this, split_mat, my_color, split_slice)
       IMPORT :: Matrix_p
       IMPLICIT NONE
       CLASS(Matrix_p), INTENT(INOUT) :: this
       CLASS(Matrix_p), INTENT(INOUT) :: split_mat
       INTEGER, INTENT(OUT) :: my_color
       LOGICAL, INTENT(OUT) :: split_slice
     END SUBROUTINE
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the actual dimension of the matrix.
  !! @param[in] this the matrix.
  !! @result dimension of the matrix;
  PURE FUNCTION GetMatrixActualDimension_p(this) RESULT(DIMENSION)
    !! Parameters
    CLASS(Matrix_p), INTENT(IN) :: this
    INTEGER :: DIMENSION
    DIMENSION = this%actual_matrix_dimension
  END FUNCTION GetMatrixActualDimension_p
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the logical dimension of the matrix.
  !! Includes padding.
  !! @param[in] this the matrix.
  !! @result dimension of the matrix;
  PURE FUNCTION GetMatrixLogicalDimension_p(this) RESULT(DIMENSION)
    !! Parameters
    CLASS(Matrix_p), INTENT(IN) :: this
    INTEGER :: DIMENSION
    DIMENSION = this%logical_matrix_dimension
  END FUNCTION GetMatrixLogicalDimension_p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE MatrixPModule
