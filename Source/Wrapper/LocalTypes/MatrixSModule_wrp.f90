!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for wrapping the sparse matrix data type.
MODULE MatrixSModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE MatrixDModule
  USE MatrixMemoryPoolModule_wrp
  USE MatrixSModule
  USE TripletListModule_wrp, ONLY : TripletList_r_wrp, TripletList_c_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int, c_char, c_bool
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for the sparse matrix data type.
  TYPE, PUBLIC :: Matrix_lsr_wrp
     TYPE(Matrix_lsr), POINTER :: DATA
  END TYPE Matrix_lsr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for the sparse matrix data type.
  TYPE, PUBLIC :: Matrix_lsc_wrp
     TYPE(Matrix_lsc), POINTER :: DATA
  END TYPE Matrix_lsc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructMatrixFromFile_lsr_wrp
  PUBLIC :: ConstructMatrixFromTripletList_lsr_wrp
  PUBLIC :: ConstructZeroMatrix_lsr_wrp
  PUBLIC :: DestructMatrix_lsr_wrp
  PUBLIC :: CopyMatrix_lsr_wrp
  PUBLIC :: GetMatrixRows_lsr_wrp
  PUBLIC :: GetMatrixColumns_lsr_wrp
  PUBLIC :: ExtractMatrixRow_lsr_wrp
  PUBLIC :: ExtractMatrixColumn_lsr_wrp
  PUBLIC :: TransposeMatrix_lsr_wrp
  PUBLIC :: PrintMatrix_lsr_wrp
  PUBLIC :: PrintMatrixF_lsr_wrp
  PUBLIC :: MatrixToTripletList_lsr_wrp
  PUBLIC :: EigenDecomposition_lsr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructMatrixFromFile_lsc_wrp
  PUBLIC :: ConstructMatrixFromTripletList_lsc_wrp
  PUBLIC :: ConstructZeroMatrix_lsc_wrp
  PUBLIC :: DestructMatrix_lsc_wrp
  PUBLIC :: CopyMatrix_lsc_wrp
  PUBLIC :: GetMatrixRows_lsc_wrp
  PUBLIC :: GetMatrixColumns_lsc_wrp
  PUBLIC :: ExtractMatrixRow_lsc_wrp
  PUBLIC :: ExtractMatrixColumn_lsc_wrp
  PUBLIC :: TransposeMatrix_lsc_wrp
  PUBLIC :: PrintMatrix_lsc_wrp
  PUBLIC :: PrintMatrixF_lsc_wrp
  PUBLIC :: MatrixToTripletList_lsc_wrp
  PUBLIC :: EigenDecomposition_lsc_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a sparse matrix by reading in a matrix market file.
  SUBROUTINE ConstructMatrixFromFile_lsr_wrp(ih_this, file_name, name_size) &
       & bind(c,name="ConstructMatrixFromFile_lsr_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    CHARACTER(kind=c_char), INTENT(IN) :: file_name(name_size)
    INTEGER(kind=c_int), INTENT(IN) :: name_size
    TYPE(Matrix_lsr_wrp) :: h_this
    !! Local Data
    CHARACTER(len=name_size) :: local_string
    INTEGER :: counter

    DO counter=1,name_size
       local_string(counter:counter) = file_name(counter)
    END DO

    ALLOCATE(h_this%data)
    CALL ConstructMatrixFromFile(h_this%data,local_string)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructMatrixFromFile_lsr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a sparse matrix from a \b SORTED triplet list.
  PURE SUBROUTINE ConstructMatrixFromTripletList_lsr_wrp(ih_this, &
       & ih_triplet_list, rows, columns) &
       & bind(c,name="ConstructMatrixFromTripletList_lsr_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_triplet_list(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: columns
    INTEGER(kind=c_int), INTENT(IN) :: rows
    !! Local Data
    TYPE(Matrix_lsr_wrp) :: h_this
    TYPE(TripletList_r_wrp)  :: h_triplet_list

    h_triplet_list = TRANSFER(ih_triplet_list,h_triplet_list)
    ALLOCATE(h_this%data)
    CALL ConstructMatrixFromTripletList(h_this%data, h_triplet_list%data, &
         & rows, columns)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructMatrixFromTripletList_lsr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a sparse matrix with zero values in it.
  PURE SUBROUTINE ConstructZeroMatrix_lsr_wrp(ih_this, rows, columns) &
       & bind(c,name="ConstructZeroMatrix_lsr_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: columns
    INTEGER(kind=c_int), INTENT(IN) :: rows
    !! Local Data
    TYPE(Matrix_lsr_wrp) :: h_this

    ALLOCATE(h_this%data)
    CALL ConstructZeroMatrix(h_this%data, rows, columns)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructZeroMatrix_lsr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Explicitly destruct a sparse matrix
  PURE SUBROUTINE DestructMatrix_lsr_wrp(ih_this) &
       & bind(c,name="DestructMatrix_lsr_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    TYPE(Matrix_lsr_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL DestructMatrix(h_this%data)
    DEALLOCATE(h_this%data)
    !ih_this = 0
  END SUBROUTINE DestructMatrix_lsr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the copy a sparse matrix routine.
  PURE SUBROUTINE CopyMatrix_lsr_wrp(ih_matA, ih_matB) &
       & bind(c,name="CopyMatrix_lsr_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matB(SIZE_wrp)
    TYPE(Matrix_lsr_wrp) :: h_matA
    TYPE(Matrix_lsr_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    CALL CopyMatrix(h_matA%data,h_matB%data)
  END SUBROUTINE CopyMatrix_lsr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the row accessor.
  PURE SUBROUTINE GetMatrixRows_lsr_wrp(ih_this, rows) &
       & bind(c,name="GetMatrixRows_lsr_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(OUT) :: rows
    TYPE(Matrix_lsr_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    rows = GetMatrixRows(h_this%data)
  END SUBROUTINE GetMatrixRows_lsr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the column accessor.
  PURE SUBROUTINE GetMatrixColumns_lsr_wrp(ih_this, columns) &
       & bind(c,name="GetMatrixColumns_lsr_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(OUT) :: columns
    TYPE(Matrix_lsr_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    columns = GetMatrixColumns(h_this%data)
  END SUBROUTINE GetMatrixColumns_lsr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract a row from the matrix into the compressed vector representation.
  !! @param[in] ih_this the matrix to extrat from.
  !! @param[in] row_number the row to extract
  !! @param[out] ih_row_out the matrix representing that row
  PURE SUBROUTINE ExtractMatrixRow_lsr_wrp(ih_this, row_number, ih_row_out) &
       & bind(c,name="ExtractMatrixRow_lsr_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: row_number
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_row_out(SIZE_wrp)
    TYPE(Matrix_lsr_wrp) :: h_this
    TYPE(Matrix_lsr_wrp) :: h_row_out

    ALLOCATE(h_row_out%data)
    h_this = TRANSFER(ih_this,h_this)
    CALL ExtractMatrixRow(h_this%data, row_number, h_row_out%data)

    ih_row_out= TRANSFER(h_row_out,ih_row_out)
  END SUBROUTINE ExtractMatrixRow_lsr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract a column from the matrix into the compressed vector representation.
  !! @param[in] ih_this the matrix to extrat from.
  !! @param[in] column_number the row to extract.
  !! @param[out] ih_column_out the matrix representing that column.
  PURE SUBROUTINE ExtractMatrixColumn_lsr_wrp(ih_this, column_number, &
       & ih_column_out) bind(c,name="ExtractMatrixColumn_lsr_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: column_number
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_column_out(SIZE_wrp)
    TYPE(Matrix_lsr_wrp) :: h_this
    TYPE(Matrix_lsr_wrp) :: h_column_out

    ALLOCATE(h_column_out%data)
    h_this = TRANSFER(ih_this,h_this)
    CALL ExtractMatrixColumn(h_this%data, column_number, h_column_out%data)

    ih_column_out= TRANSFER(h_column_out, ih_column_out)
  END SUBROUTINE ExtractMatrixColumn_lsr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the matrix transpose function.
  PURE SUBROUTINE TransposeMatrix_lsr_wrp(ih_matA, ih_matAT) &
       & bind(c,name="TransposeMatrix_lsr_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matAT(SIZE_wrp)
    TYPE(Matrix_lsr_wrp) :: h_matA
    TYPE(Matrix_lsr_wrp) :: h_matAT

    h_matA  = TRANSFER(ih_matA,h_matA)
    h_matAT = TRANSFER(ih_matAT,h_matAT)
    CALL TransposeMatrix(h_matA%data,h_matAT%data)
  END SUBROUTINE TransposeMatrix_lsr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Warp the routine that prints out a sparse matrix to file.
  SUBROUTINE PrintMatrixF_lsr_wrp(ih_this, file_name, name_size) &
       & bind(c,name="PrintMatrixF_lsr_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    CHARACTER(kind=c_char), INTENT(IN) :: file_name(name_size)
    INTEGER(kind=c_int), INTENT(IN) :: name_size
    TYPE(Matrix_lsr_wrp) :: h_this
    !! Local Data
    CHARACTER(len=name_size) :: local_string
    INTEGER :: counter

    DO counter=1,name_size
       local_string(counter:counter) = file_name(counter)
    END DO

    h_this = TRANSFER(ih_this,h_this)
    CALL PrintMatrix(h_this%data,local_string)
  END SUBROUTINE PrintMatrixF_lsr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Warp the routine that prints the sparse matrix to the console.
  SUBROUTINE PrintMatrix_lsr_wrp(ih_this) bind(c,name="PrintMatrix_lsr_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    TYPE(Matrix_lsr_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL PrintMatrix(h_this%data)
  END SUBROUTINE PrintMatrix_lsr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the routine that constructs a triplet list from a matrix.
  SUBROUTINE MatrixToTripletList_lsr_wrp(ih_this, ih_triplet_list) &
       & bind(c,name="MatrixToTripletList_lsr_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(OUT)   :: ih_triplet_list(SIZE_wrp)
    TYPE(Matrix_lsr_wrp) :: h_this
    TYPE(TripletList_r_wrp)  :: h_triplet_list

    h_this = TRANSFER(ih_this,h_this)
    ALLOCATE(h_triplet_list%data)

    CALL MatrixToTripletList(h_this%data,h_triplet_list%data)

    ih_triplet_list = TRANSFER(ih_triplet_list,ih_triplet_list)
  END SUBROUTINE MatrixToTripletList_lsr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the dense eigendecompsition routine.
  SUBROUTINE EigenDecomposition_lsr_wrp(ih_matA, ih_matV, threshold) &
       & bind(c,name="EigenDecomposition_lsr_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matV(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: threshold
    TYPE(Matrix_lsr_wrp) :: h_matA
    TYPE(Matrix_lsr_wrp) :: h_matV
    !! Temporary Dense Matrices
    TYPE(Matrix_ldr) :: dense_a, dense_v

    h_matA  = TRANSFER(ih_matA,h_matA)
    h_matV = TRANSFER(ih_matV,h_matV)

    CALL ConstructMatrixDFromS(h_matA%data, dense_a)
    CALL EigenDecomposition(dense_a, dense_v)
    CALL ConstructMatrixSFromD(dense_v,h_matV%data,threshold)

    !! Cleanup
    CALL DestructMatrix(dense_a)
    CALL DestructMatrix(dense_v)
  END SUBROUTINE EigenDecomposition_lsr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a sparse matrix by reading in a matrix market file.
  SUBROUTINE ConstructMatrixFromFile_lsc_wrp(ih_this, file_name, name_size) &
       & bind(c,name="ConstructMatrixFromFile_lsc_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    CHARACTER(kind=c_char), INTENT(IN) :: file_name(name_size)
    INTEGER(kind=c_int), INTENT(IN) :: name_size
    TYPE(Matrix_lsc_wrp) :: h_this
    !! Local Data
    CHARACTER(len=name_size) :: local_string
    INTEGER :: counter

    DO counter=1,name_size
       local_string(counter:counter) = file_name(counter)
    END DO

    ALLOCATE(h_this%data)
    CALL ConstructMatrixFromFile(h_this%data,local_string)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructMatrixFromFile_lsc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a sparse matrix from a \b SORTED triplet list.
  PURE SUBROUTINE ConstructMatrixFromTripletList_lsc_wrp(ih_this, &
       & ih_triplet_list, rows, columns) &
       & bind(c,name="ConstructMatrixFromTripletList_lsc_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_triplet_list(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: columns
    INTEGER(kind=c_int), INTENT(IN) :: rows
    !! Local Data
    TYPE(Matrix_lsc_wrp) :: h_this
    TYPE(TripletList_c_wrp)  :: h_triplet_list

    h_triplet_list = TRANSFER(ih_triplet_list,h_triplet_list)
    ALLOCATE(h_this%data)
    CALL ConstructMatrixFromTripletList(h_this%data, h_triplet_list%data, &
         & rows, columns)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructMatrixFromTripletList_lsc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a sparse matrix with zero values in it.
  PURE SUBROUTINE ConstructZeroMatrix_lsc_wrp(ih_this, rows, columns) &
       & bind(c,name="ConstructZeroMatrix_lsc_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: columns
    INTEGER(kind=c_int), INTENT(IN) :: rows
    !! Local Data
    TYPE(Matrix_lsc_wrp) :: h_this

    ALLOCATE(h_this%data)
    CALL ConstructZeroMatrix(h_this%data, rows, columns)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructZeroMatrix_lsc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Explicitly destruct a sparse matrix
  PURE SUBROUTINE DestructMatrix_lsc_wrp(ih_this) &
       & bind(c,name="DestructMatrix_lsc_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    TYPE(Matrix_lsc_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL DestructMatrix(h_this%data)
    DEALLOCATE(h_this%data)
    !ih_this = 0
  END SUBROUTINE DestructMatrix_lsc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the copy a sparse matrix routine.
  PURE SUBROUTINE CopyMatrix_lsc_wrp(ih_matA, ih_matB) &
       & bind(c,name="CopyMatrix_lsc_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matB(SIZE_wrp)
    TYPE(Matrix_lsc_wrp) :: h_matA
    TYPE(Matrix_lsc_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    CALL CopyMatrix(h_matA%data,h_matB%data)
  END SUBROUTINE CopyMatrix_lsc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the row accessor.
  PURE SUBROUTINE GetMatrixRows_lsc_wrp(ih_this, rows) &
       & bind(c,name="GetMatrixRows_lsc_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(OUT) :: rows
    TYPE(Matrix_lsc_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    rows = GetMatrixRows(h_this%data)
  END SUBROUTINE GetMatrixRows_lsc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the column accessor.
  PURE SUBROUTINE GetMatrixColumns_lsc_wrp(ih_this, columns) &
       & bind(c,name="GetMatrixColumns_lsc_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(OUT) :: columns
    TYPE(Matrix_lsc_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    columns = GetMatrixColumns(h_this%data)
  END SUBROUTINE GetMatrixColumns_lsc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract a row from the matrix into the compressed vector representation.
  !! @param[in] ih_this the matrix to extrat from.
  !! @param[in] row_number the row to extract
  !! @param[out] ih_row_out the matrix representing that row
  PURE SUBROUTINE ExtractMatrixRow_lsc_wrp(ih_this, row_number, ih_row_out) &
       & bind(c,name="ExtractMatrixRow_lsc_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: row_number
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_row_out(SIZE_wrp)
    TYPE(Matrix_lsc_wrp) :: h_this
    TYPE(Matrix_lsc_wrp) :: h_row_out

    ALLOCATE(h_row_out%data)
    h_this = TRANSFER(ih_this,h_this)
    CALL ExtractMatrixRow(h_this%data, row_number, h_row_out%data)

    ih_row_out= TRANSFER(h_row_out,ih_row_out)
  END SUBROUTINE ExtractMatrixRow_lsc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract a column from the matrix into the compressed vector representation.
  !! @param[in] ih_this the matrix to extrat from.
  !! @param[in] column_number the row to extract.
  !! @param[out] ih_column_out the matrix representing that column.
  PURE SUBROUTINE ExtractMatrixColumn_lsc_wrp(ih_this, column_number, &
       & ih_column_out) bind(c,name="ExtractMatrixColumn_lsc_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: column_number
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_column_out(SIZE_wrp)
    TYPE(Matrix_lsc_wrp) :: h_this
    TYPE(Matrix_lsc_wrp) :: h_column_out

    ALLOCATE(h_column_out%data)
    h_this = TRANSFER(ih_this,h_this)
    CALL ExtractMatrixColumn(h_this%data, column_number, h_column_out%data)

    ih_column_out= TRANSFER(h_column_out, ih_column_out)
  END SUBROUTINE ExtractMatrixColumn_lsc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the matrix transpose function.
  PURE SUBROUTINE TransposeMatrix_lsc_wrp(ih_matA, ih_matAT) &
       & bind(c,name="TransposeMatrix_lsc_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matAT(SIZE_wrp)
    TYPE(Matrix_lsc_wrp) :: h_matA
    TYPE(Matrix_lsc_wrp) :: h_matAT

    h_matA  = TRANSFER(ih_matA,h_matA)
    h_matAT = TRANSFER(ih_matAT,h_matAT)
    CALL TransposeMatrix(h_matA%data,h_matAT%data)
  END SUBROUTINE TransposeMatrix_lsc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Warp the routine that prints out a sparse matrix to file.
  SUBROUTINE PrintMatrixF_lsc_wrp(ih_this, file_name, name_size) &
       & bind(c,name="PrintMatrixF_lsc_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    CHARACTER(kind=c_char), INTENT(IN) :: file_name(name_size)
    INTEGER(kind=c_int), INTENT(IN) :: name_size
    TYPE(Matrix_lsc_wrp) :: h_this
    !! Local Data
    CHARACTER(len=name_size) :: local_string
    INTEGER :: counter

    DO counter=1,name_size
       local_string(counter:counter) = file_name(counter)
    END DO

    h_this = TRANSFER(ih_this,h_this)
    CALL PrintMatrix(h_this%data,local_string)
  END SUBROUTINE PrintMatrixF_lsc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Warp the routine that prints the sparse matrix to the console.
  SUBROUTINE PrintMatrix_lsc_wrp(ih_this) bind(c,name="PrintMatrix_lsc_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    TYPE(Matrix_lsc_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL PrintMatrix(h_this%data)
  END SUBROUTINE PrintMatrix_lsc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the routine that constructs a triplet list from a matrix.
  SUBROUTINE MatrixToTripletList_lsc_wrp(ih_this, ih_triplet_list) &
       & bind(c,name="MatrixToTripletList_lsc_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(OUT)   :: ih_triplet_list(SIZE_wrp)
    TYPE(Matrix_lsc_wrp) :: h_this
    TYPE(TripletList_c_wrp)  :: h_triplet_list

    h_this = TRANSFER(ih_this,h_this)
    ALLOCATE(h_triplet_list%data)

    CALL MatrixToTripletList(h_this%data,h_triplet_list%data)

    ih_triplet_list = TRANSFER(ih_triplet_list,ih_triplet_list)
  END SUBROUTINE MatrixToTripletList_lsc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the dense eigendecompsition routine.
  SUBROUTINE EigenDecomposition_lsc_wrp(ih_matA, ih_matV, threshold) &
       & bind(c,name="EigenDecomposition_lsc_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matV(SIZE_wrp)
    REAL(NTREAL), INTENT(IN) :: threshold
    TYPE(Matrix_lsc_wrp) :: h_matA
    TYPE(Matrix_lsc_wrp) :: h_matV
    !! Temporary Dense Matrices
    TYPE(Matrix_ldc) :: dense_a, dense_v

    h_matA  = TRANSFER(ih_matA,h_matA)
    h_matV = TRANSFER(ih_matV,h_matV)

    CALL ConstructMatrixDFromS(h_matA%data, dense_a)
    CALL EigenDecomposition(dense_a, dense_v)
    CALL ConstructMatrixSFromD(dense_v,h_matV%data,threshold)

    !! Cleanup
    CALL DestructMatrix(dense_a)
    CALL DestructMatrix(dense_v)
  END SUBROUTINE EigenDecomposition_lsc_wrp
END MODULE MatrixSModule_wrp
