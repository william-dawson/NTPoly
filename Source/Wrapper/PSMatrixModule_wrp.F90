!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for wrapping a Distributed Sparse Matrix.
MODULE PSMatrixModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE PermutationModule_wrp, ONLY : Permutation_wrp
  USE ProcessGridModule_wrp, ONLY : ProcessGrid_wrp
  USE PSMatrixModule
  USE TripletListModule_wrp, ONLY : TripletList_r_wrp, TripletList_c_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int, c_char, c_bool
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for the distributed sparse matrix data type.
  TYPE, PUBLIC :: Matrix_ps_wrp
     TYPE(Matrix_ps), POINTER :: DATA
  END TYPE Matrix_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructEmptyMatrix_ps_wrp
  PUBLIC :: CopyMatrix_ps_wrp
  PUBLIC :: DestructMatrix_ps_wrp
  PUBLIC :: ConstructMatrixFromMatrixMarket_ps_wrp
  PUBLIC :: ConstructMatrixFromBinary_ps_wrp
  PUBLIC :: WriteMatrixToBinary_ps_wrp
  PUBLIC :: WriteMatrixToMatrixMarket_ps_wrp
  PUBLIC :: FillMatrixFromTripletList_psr_wrp
  PUBLIC :: FillMatrixFromTripletList_psc_wrp
  PUBLIC :: FillMatrixPermutation_ps_wrp
  PUBLIC :: FillMatrixIdentity_ps_wrp
  PUBLIC :: GetMatrixActualDimension_ps_wrp
  PUBLIC :: GetMatrixLogicalDimension_ps_wrp
  PUBLIC :: GetMatrixTripletList_psr_wrp
  PUBLIC :: GetMatrixTripletList_psc_wrp
  PUBLIC :: GetMatrixBlock_psr_wrp
  PUBLIC :: GetMatrixBlock_psc_wrp
  PUBLIC :: GetMatrixSlice_wrp
  PUBLIC :: TransposeMatrix_ps_wrp
  PUBLIC :: ConjugateMatrix_ps_wrp
  PUBLIC :: ResizeMatrix_ps_wrp
  PUBLIC :: GetMatrixProcessGrid_ps_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the constructor of an empty sparse, distributed, matrix.
  SUBROUTINE ConstructEmptyMatrix_ps_wrp(ih_this,matrix_dim) &
       & BIND(c,NAME="ConstructEmptyMatrix_ps_wrp")
    INTEGER(KIND=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(KIND=c_int), INTENT(IN) :: matrix_dim
    TYPE(Matrix_ps_wrp) :: h_this

    ALLOCATE(h_this%data)
    CALL ConstructEmptyMatrix(h_this%data,matrix_dim)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructEmptyMatrix_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the constructor of an empty sparse, distributed, matrix.
  SUBROUTINE ConstructEmptyMatrixPG_ps_wrp(ih_this,matrix_dim,ih_grid) &
       & BIND(c,NAME="ConstructEmptyMatrixPG_ps_wrp")
    INTEGER(KIND=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(KIND=c_int), INTENT(IN) :: matrix_dim
    INTEGER(KIND=c_int), INTENT(IN) :: ih_grid(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_this
    TYPE(ProcessGrid_wrp) :: h_grid

    h_grid = TRANSFER(ih_grid,h_grid)
    ALLOCATE(h_this%data)
    CALL ConstructEmptyMatrix(h_this%data,matrix_dim,h_grid%data)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructEmptyMatrixPG_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct distributed sparse matrix from a matrix market file in parallel.
  SUBROUTINE ConstructMatrixFromMatrixMarket_ps_wrp(ih_this, file_name, &
       & name_size) BIND(c,NAME="ConstructMatrixFromMatrixMarket_ps_wrp")
    INTEGER(KIND=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    CHARACTER(KIND=c_char), INTENT(IN) :: file_name(name_size)
    INTEGER(KIND=c_int), INTENT(IN) :: name_size
    TYPE(Matrix_ps_wrp) :: h_this
    !! Local Data
    CHARACTER(len=name_size) :: local_string
    INTEGER :: counter

    DO counter=1,name_size
       local_string(counter:counter) = file_name(counter)
    END DO

    ALLOCATE(h_this%data)
    CALL ConstructMatrixFromMatrixMarket(h_this%data,local_string)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructMatrixFromMatrixMarket_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct distributed sparse matrix from a matrix market file in parallel.
  SUBROUTINE ConstructMatrixFromMatrixMarketPG_ps_wrp(ih_this, file_name, &
       & name_size, ih_grid) &
       & BIND(c,NAME="ConstructMatrixFromMatrixMarketPG_ps_wrp")
    INTEGER(KIND=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    CHARACTER(KIND=c_char), INTENT(IN) :: file_name(name_size)
    INTEGER(KIND=c_int), INTENT(IN) :: name_size
    INTEGER(KIND=c_int), INTENT(IN) :: ih_grid(SIZE_wrp)
    !! Local Data
    TYPE(Matrix_ps_wrp) :: h_this
    TYPE(ProcessGrid_wrp) :: h_grid
    CHARACTER(len=name_size) :: local_string
    INTEGER :: counter

    DO counter=1,name_size
       local_string(counter:counter) = file_name(counter)
    END DO

    ALLOCATE(h_this%data)
    h_grid = TRANSFER(ih_grid,h_grid)
    CALL ConstructMatrixFromMatrixMarket(h_this%data,local_string,h_grid%data)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructMatrixFromMatrixMarketPG_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a distributed sparse matrix from a binary file in parallel.
  SUBROUTINE ConstructMatrixFromBinary_ps_wrp(ih_this,file_name,name_size) &
       & BIND(c,NAME="ConstructMatrixFromBinary_ps_wrp")
    INTEGER(KIND=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    CHARACTER(KIND=c_char), INTENT(IN) :: file_name(name_size)
    INTEGER(KIND=c_int), INTENT(IN) :: name_size
    !! Local Data
    TYPE(Matrix_ps_wrp) :: h_this
    CHARACTER(len=name_size) :: local_string
    INTEGER :: counter

    DO counter=1,name_size
       local_string(counter:counter) = file_name(counter)
    END DO

    ALLOCATE(h_this%data)
    CALL ConstructMatrixFromBinary(h_this%data,local_string)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructMatrixFromBinary_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a distributed sparse matrix from a binary file in parallel.
  SUBROUTINE ConstructMatrixFromBinaryPG_ps_wrp(ih_this, file_name, &
       & name_size, ih_grid) &
       & BIND(c,NAME="ConstructMatrixFromBinaryPG_ps_wrp")
    INTEGER(KIND=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    CHARACTER(KIND=c_char), INTENT(IN) :: file_name(name_size)
    INTEGER(KIND=c_int), INTENT(IN) :: name_size
    INTEGER(KIND=c_int), INTENT(IN) :: ih_grid(SIZE_wrp)
    !! Local Data
    TYPE(Matrix_ps_wrp) :: h_this
    TYPE(ProcessGrid_wrp) :: h_grid
    CHARACTER(len=name_size) :: local_string
    INTEGER :: counter

    DO counter=1,name_size
       local_string(counter:counter) = file_name(counter)
    END DO

    ALLOCATE(h_this%data)
    h_grid = TRANSFER(ih_grid,h_grid)
    CALL ConstructMatrixFromBinary(h_this%data,local_string,h_grid%data)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructMatrixFromBinaryPG_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy a distributed sparse matrix in a safe way.
  SUBROUTINE CopyMatrix_ps_wrp(ih_matA,ih_matB) &
       & BIND(c,NAME="CopyMatrix_ps_wrp")
    INTEGER(KIND=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(KIND=c_int), INTENT(INOUT) :: ih_matB(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_matA
    TYPE(Matrix_ps_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    CALL CopyMatrix(h_matA%data,h_matB%data)
  END SUBROUTINE CopyMatrix_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a distributed sparse matrix
  PURE SUBROUTINE DestructMatrix_ps_wrp(ih_this) &
       & BIND(c,NAME="DestructMatrix_ps_wrp")
    INTEGER(KIND=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL DestructMatrix(h_this%data)
    DEALLOCATE(h_this%data)
  END SUBROUTINE DestructMatrix_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Save a distributed sparse matrix to a file.
  SUBROUTINE WriteMatrixToBinary_ps_wrp(ih_this,file_name,name_size) &
       & BIND(c,NAME="WriteMatrixToBinary_ps_wrp")
    INTEGER(KIND=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    CHARACTER(KIND=c_char), INTENT(IN) :: file_name(name_size)
    INTEGER(KIND=c_int), INTENT(IN) :: name_size
    TYPE(Matrix_ps_wrp) :: h_this
    !! Local Data
    CHARACTER(len=name_size) :: local_string
    INTEGER :: counter

    DO counter=1,name_size
       local_string(counter:counter) = file_name(counter)
    END DO

    h_this = TRANSFER(ih_this,h_this)
    CALL WriteMatrixToBinary(h_this%data,local_string)
  END SUBROUTINE WriteMatrixToBinary_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Save a distributed sparse matrix to a matrix market file.
  SUBROUTINE WriteMatrixToMatrixMarket_ps_wrp(ih_this,file_name,name_size) &
       & BIND(c,NAME="WriteMatrixToMatrixMarket_ps_wrp")
    INTEGER(KIND=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    CHARACTER(KIND=c_char), INTENT(IN) :: file_name(name_size)
    INTEGER(KIND=c_int), INTENT(IN) :: name_size
    TYPE(Matrix_ps_wrp) :: h_this
    !! Local Data
    CHARACTER(len=name_size) :: local_string
    INTEGER :: counter

    DO counter=1,name_size
       local_string(counter:counter) = file_name(counter)
    END DO

    h_this = TRANSFER(ih_this,h_this)
    CALL WriteMatrixToMatrixMarket(h_this%data,local_string)
  END SUBROUTINE WriteMatrixToMatrixMarket_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This routine fills in a matrix based on local triplet lists.
  SUBROUTINE FillMatrixFromTripletList_psr_wrp(ih_this, ih_triplet_list) &
       & BIND(c,NAME="FillMatrixFromTripletList_psr_wrp")
    INTEGER(KIND=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(KIND=c_int), INTENT(IN) :: ih_triplet_list(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_this
    TYPE(TripletList_r_wrp) :: h_triplet_list

    h_this = TRANSFER(ih_this,h_this)
    h_triplet_list = TRANSFER(ih_triplet_list,h_triplet_list)
    CALL FillMatrixFromTripletList(h_this%data, h_triplet_list%data)
  END SUBROUTINE FillMatrixFromTripletList_psr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This routine fills in a matrix based on local triplet lists.
  SUBROUTINE FillMatrixFromTripletList_psc_wrp(ih_this, ih_triplet_list) &
       & BIND(c,NAME="FillMatrixFromTripletList_psc_wrp")
    INTEGER(KIND=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(KIND=c_int), INTENT(IN) :: ih_triplet_list(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_this
    TYPE(TripletList_c_wrp) :: h_triplet_list

    h_this = TRANSFER(ih_this,h_this)
    h_triplet_list = TRANSFER(ih_triplet_list,h_triplet_list)
    CALL FillMatrixFromTripletList(h_this%data, h_triplet_list%data)
  END SUBROUTINE FillMatrixFromTripletList_psc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Fill in the values of a distributed matrix with the identity matrix.
  SUBROUTINE FillMatrixIdentity_ps_wrp(ih_this) &
       & BIND(c,NAME="FillMatrixIdentity_ps_wrp")
    INTEGER(KIND=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL FillMatrixIdentity(h_this%data)
  END SUBROUTINE FillMatrixIdentity_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Fill in the values of a distributed matrix with a permutation matrix.
  SUBROUTINE FillMatrixPermutation_ps_wrp(ih_this, ih_permutation, &
       & permute_rows) BIND(c,NAME="FillMatrixPermutation_ps_wrp")
    INTEGER(KIND=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(KIND=c_int), INTENT(IN) :: ih_permutation(SIZE_wrp)
    LOGICAL(KIND=c_bool), INTENT(IN) :: permute_rows
    TYPE(Matrix_ps_wrp) :: h_this
    TYPE(Permutation_wrp) :: h_permutation

    h_this = TRANSFER(ih_this,h_this)
    h_permutation = TRANSFER(ih_permutation,h_permutation)

    CALL FillMatrixPermutation(h_this%data,h_permutation%data%index_lookup, &
         & permute_rows_in = LOGICAL(permute_rows))
  END SUBROUTINE FillMatrixPermutation_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the Actual Dimension accessor.
  PURE SUBROUTINE GetMatrixActualDimension_ps_wrp(ih_this, mat_dimension) &
       & BIND(c,NAME="GetMatrixActualDimension_ps_wrp")
    INTEGER(KIND=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(KIND=c_int), INTENT(OUT) :: mat_dimension
    TYPE(Matrix_ps_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    mat_dimension = GetMatrixActualDimension(h_this%data)
  END SUBROUTINE GetMatrixActualDimension_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the Logical Dimension accessor.
  PURE SUBROUTINE GetMatrixLogicalDimension_ps_wrp(ih_this, mat_dimension) &
       & BIND(c,NAME="GetMatrixLogicalDimension_ps_wrp")
    INTEGER(KIND=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(KIND=c_int), INTENT(OUT) :: mat_dimension
    TYPE(Matrix_ps_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    mat_dimension = GetMatrixLogicalDimension(h_this%data)
  END SUBROUTINE GetMatrixLogicalDimension_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extracts a triplet list of the data that is stored on this process.
  SUBROUTINE GetMatrixTripletList_psr_wrp(ih_this, ih_triplet_list) &
       & BIND(c,NAME="GetMatrixTripletList_psr_wrp")
    INTEGER(KIND=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(KIND=c_int), INTENT(INOUT) :: ih_triplet_list(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_this
    TYPE(TripletList_r_wrp) :: h_triplet_list

    h_this = TRANSFER(ih_this,h_this)
    h_triplet_list = TRANSFER(ih_triplet_list,h_triplet_list)
    CALL GetMatrixTripletList(h_this%data,h_triplet_list%data)
  END SUBROUTINE GetMatrixTripletList_psr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extracts a triplet list of the data that is stored on this process.
  SUBROUTINE GetMatrixTripletList_psc_wrp(ih_this, ih_triplet_list) &
       & BIND(c,NAME="GetMatrixTripletList_psc_wrp")
    INTEGER(KIND=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(KIND=c_int), INTENT(INOUT) :: ih_triplet_list(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_this
    TYPE(TripletList_c_wrp) :: h_triplet_list

    h_this = TRANSFER(ih_this,h_this)
    h_triplet_list = TRANSFER(ih_triplet_list,h_triplet_list)
    CALL GetMatrixTripletList(h_this%data,h_triplet_list%data)
  END SUBROUTINE GetMatrixTripletList_psc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract an arbitrary block of a matrix into a triplet list.
  SUBROUTINE GetMatrixBlock_psr_wrp(ih_this, ih_triplet_list, start_row, &
       & end_row, start_column, end_column) BIND(c,NAME="GetMatrixBlock_psr_wrp")
    INTEGER(KIND=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(KIND=c_int), INTENT(INOUT) :: ih_triplet_list(SIZE_wrp)
    INTEGER(KIND=c_int), INTENT(IN) :: start_row, end_row
    INTEGER(KIND=c_int), INTENT(IN) :: start_column, end_column
    TYPE(Matrix_ps_wrp) :: h_this
    TYPE(TripletList_r_wrp) :: h_triplet_list

    h_this = TRANSFER(ih_this,h_this)
    h_triplet_list = TRANSFER(ih_triplet_list,h_triplet_list)
    CALL GetMatrixBlock(h_this%data,h_triplet_list%data, start_row, end_row,&
         & start_column, end_column)
  END SUBROUTINE GetMatrixBlock_psr_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract an arbitrary block of a matrix into a triplet list.
  SUBROUTINE GetMatrixBlock_psc_wrp(ih_this, ih_triplet_list, start_row, &
       & end_row, start_column, end_column) BIND(c,NAME="GetMatrixBlock_psc_wrp")
    INTEGER(KIND=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(KIND=c_int), INTENT(INOUT) :: ih_triplet_list(SIZE_wrp)
    INTEGER(KIND=c_int), INTENT(IN) :: start_row, end_row
    INTEGER(KIND=c_int), INTENT(IN) :: start_column, end_column
    TYPE(Matrix_ps_wrp) :: h_this
    TYPE(TripletList_c_wrp) :: h_triplet_list

    h_this = TRANSFER(ih_this,h_this)
    h_triplet_list = TRANSFER(ih_triplet_list,h_triplet_list)
    CALL GetMatrixBlock(h_this%data,h_triplet_list%data, start_row, end_row,&
         & start_column, end_column)
  END SUBROUTINE GetMatrixBlock_psc_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy an arbitrary slice from a matrix into a new smaller matrix.
  SUBROUTINE GetMatrixSlice_wrp(ih_this, ih_submatrix, start_row, &
       & end_row, start_column, end_column) BIND(c,NAME="GetMatrixSlice_wrp")
    INTEGER(KIND=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(KIND=c_int), INTENT(INOUT) :: ih_submatrix(SIZE_wrp)
    INTEGER(KIND=c_int), INTENT(IN) :: start_row, end_row
    INTEGER(KIND=c_int), INTENT(IN) :: start_column, end_column
    TYPE(Matrix_ps_wrp) :: h_this
    TYPE(Matrix_ps_wrp) :: h_submatrix

    h_this = TRANSFER(ih_this,h_this)
    h_submatrix = TRANSFER(ih_submatrix,h_submatrix)
    CALL GetMatrixSlice(h_this%data,h_submatrix%data, start_row, end_row,&
         & start_column, end_column)
  END SUBROUTINE GetMatrixSlice_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Transpose a sparse matrix.
  SUBROUTINE TransposeMatrix_ps_wrp(ih_matA,ih_transmat) &
       & BIND(c,NAME="TransposeMatrix_ps_wrp")
    INTEGER(KIND=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(KIND=c_int), INTENT(INOUT) :: ih_transmat(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_matA
    TYPE(Matrix_ps_wrp) :: h_transmat

    h_matA = TRANSFER(ih_matA,h_matA)
    h_transmat = TRANSFER(ih_transmat,h_transmat)
    CALL TransposeMatrix(h_matA%data,h_transmat%data)
  END SUBROUTINE TransposeMatrix_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the matrix conjugate function.
  PURE SUBROUTINE ConjugateMatrix_ps_wrp(ih_matA) &
       & BIND(c,NAME="ConjugateMatrix_ps_wrp")
    INTEGER(KIND=c_int), INTENT(INOUT) :: ih_matA(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_matA

    h_matA  = TRANSFER(ih_matA,h_matA)
    CALL ConjugateMatrix(h_matA%data)
  END SUBROUTINE ConjugateMatrix_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the matrix resize function.
  SUBROUTINE ResizeMatrix_ps_wrp(ih_this, new_size) &
       & BIND(c,NAME="ResizeMatrix_ps_wrp")
    INTEGER(KIND=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(KIND=c_int), INTENT(IN) :: new_size
    TYPE(Matrix_ps_wrp) :: h_this

    h_this  = TRANSFER(ih_this,h_this)
    CALL ResizeMatrix(h_this%data, new_size)
  END SUBROUTINE ResizeMatrix_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Returns a handle to the process grid this matrix is distributed on.
  SUBROUTINE GetMatrixProcessGrid_ps_wrp(ih_this, ih_grid) &
       & BIND(c,NAME="GetMatrixProcessGrid_ps_wrp")
    INTEGER(KIND=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(KIND=c_int), INTENT(INOUT) :: ih_grid(SIZE_wrp)
    TYPE(Matrix_ps_wrp) :: h_this
    TYPE(ProcessGrid_wrp) :: h_grid

    h_this = TRANSFER(ih_this,h_this)
    h_grid%data => h_this%data%process_grid
    ih_grid = TRANSFER(h_grid, ih_grid)
  END SUBROUTINE GetMatrixProcessGrid_ps_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE PSMatrixModule_wrp
