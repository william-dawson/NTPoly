!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for wrapping a Distributed Sparse Matrix.
MODULE DistributedSparseMatrixModule_wrp
  USE DataTypesModule, ONLY : NTREAL
  USE DistributedSparseMatrixAlgebraModule
  USE DistributedSparseMatrixModule
  USE PermutationModule_wrp, ONLY : Permutation_wrp
  USE TripletListModule_wrp, ONLY : TripletList_wrp
  USE WrapperModule, ONLY : SIZE_wrp
  USE ISO_C_BINDING, ONLY : c_int, c_char, c_bool
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for the distributed sparse matrix data type.
  TYPE, PUBLIC :: DistributedSparseMatrix_wrp
     TYPE(DistributedSparseMatrix_t), POINTER :: DATA
  END TYPE DistributedSparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructEmptyDistributedSparseMatrix_wrp
  PUBLIC :: CopyDistributedSparseMatrix_wrp
  PUBLIC :: DestructDistributedSparseMatrix_wrp
  PUBLIC :: ConstructFromMatrixMarket_wrp
  PUBLIC :: ConstructFromBinary_wrp
  PUBLIC :: WriteToBinary_wrp
  PUBLIC :: WriteToMatrixMarket_wrp
  PUBLIC :: FillFromTripletList_wrp
  PUBLIC :: FillDistributedPermutation_wrp
  PUBLIC :: FillDistributedIdentity_wrp
  PUBLIC :: GetActualDimension_wrp
  PUBLIC :: GetLogicalDimension_wrp
  PUBLIC :: GetTripletList_wrp
  PUBLIC :: GetMatrixBlock_wrp
  PUBLIC :: TransposeDistributedSparseMatrix_wrp
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the constructor of an empty sparse, distributed, matrix.
  PURE SUBROUTINE ConstructEmptyDistributedSparseMatrix_wrp(ih_this,matrix_dim) &
       & bind(c,name="ConstructEmptyDistributedSparseMatrix_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: matrix_dim
    TYPE(DistributedSparseMatrix_wrp) :: h_this

    ALLOCATE(h_this%data)
    CALL ConstructEmptyDistributedSparseMatrix(h_this%data,matrix_dim)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructEmptyDistributedSparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct distributed sparse matrix from a matrix market file in parallel.
  SUBROUTINE ConstructFromMatrixMarket_wrp(ih_this,file_name,name_size) &
       & bind(c,name="ConstructFromMatrixMarket_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    CHARACTER(kind=c_char), INTENT(IN) :: file_name(name_size)
    INTEGER(kind=c_int), INTENT(IN) :: name_size
    TYPE(DistributedSparseMatrix_wrp) :: h_this
    !! Local Data
    CHARACTER(len=name_size) :: local_string
    INTEGER :: counter

    DO counter=1,name_size
       local_string(counter:counter) = file_name(counter)
    END DO

    ALLOCATE(h_this%data)
    CALL ConstructFromMatrixMarket(h_this%data,local_string)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructFromMatrixMarket_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a distributed sparse matrix from a binary file in parallel.
  SUBROUTINE ConstructFromBinary_wrp(ih_this,file_name,name_size) &
       & bind(c,name="ConstructFromBinary_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    CHARACTER(kind=c_char), INTENT(IN) :: file_name(name_size)
    INTEGER(kind=c_int), INTENT(IN) :: name_size
    TYPE(DistributedSparseMatrix_wrp) :: h_this
    !! Local Data
    CHARACTER(len=name_size) :: local_string
    INTEGER :: counter

    DO counter=1,name_size
       local_string(counter:counter) = file_name(counter)
    END DO

    ALLOCATE(h_this%data)
    CALL ConstructFromBinary(h_this%data,local_string)
    ih_this = TRANSFER(h_this,ih_this)
  END SUBROUTINE ConstructFromBinary_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy a distributed sparse matrix in a safe way.
  PURE SUBROUTINE CopyDistributedSparseMatrix_wrp(ih_matA,ih_matB) &
       & bind(c,name="CopyDistributedSparseMatrix_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_matB(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_matA
    TYPE(DistributedSparseMatrix_wrp) :: h_matB

    h_matA = TRANSFER(ih_matA,h_matA)
    h_matB = TRANSFER(ih_matB,h_matB)
    CALL CopyDistributedSparseMatrix(h_matA%data,h_matB%data)
  END SUBROUTINE CopyDistributedSparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a distributed sparse matrix
  PURE SUBROUTINE DestructDistributedSparseMatrix_wrp(ih_this) &
       & bind(c,name="DestructDistributedSparseMatrix_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL DestructDistributedSparseMatrix(h_this%data)
    DEALLOCATE(h_this%data)
    !ih_this = 0
  END SUBROUTINE DestructDistributedSparseMatrix_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Save a distributed sparse matrix to a file.
  SUBROUTINE WriteToBinary_wrp(ih_this,file_name,name_size) &
       & bind(c,name="WriteToBinary_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    CHARACTER(kind=c_char), INTENT(IN) :: file_name(name_size)
    INTEGER(kind=c_int), INTENT(IN) :: name_size
    TYPE(DistributedSparseMatrix_wrp) :: h_this
    !! Local Data
    CHARACTER(len=name_size) :: local_string
    INTEGER :: counter

    DO counter=1,name_size
       local_string(counter:counter) = file_name(counter)
    END DO

    h_this = TRANSFER(ih_this,h_this)
    CALL WriteToBinary(h_this%data,local_string)
  END SUBROUTINE WriteToBinary_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Save a distributed sparse matrix to a matrix market file.
  SUBROUTINE WriteToMatrixMarket_wrp(ih_this,file_name,name_size) &
       & bind(c,name="WriteToMatrixMarket_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    CHARACTER(kind=c_char), INTENT(IN) :: file_name(name_size)
    INTEGER(kind=c_int), INTENT(IN) :: name_size
    TYPE(DistributedSparseMatrix_wrp) :: h_this
    !! Local Data
    CHARACTER(len=name_size) :: local_string
    INTEGER :: counter

    DO counter=1,name_size
       local_string(counter:counter) = file_name(counter)
    END DO

    h_this = TRANSFER(ih_this,h_this)
    CALL WriteToMatrixMarket(h_this%data,local_string)
  END SUBROUTINE WriteToMatrixMarket_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This routine fills in a matrix based on local triplet lists.
  SUBROUTINE FillFromTripletList_wrp(ih_this, ih_triplet_list) &
       & bind(c,name="FillFromTripletList_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_triplet_list(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_this
    TYPE(TripletList_wrp) :: h_triplet_list

    h_this = TRANSFER(ih_this,h_this)
    h_triplet_list = TRANSFER(ih_triplet_list,h_triplet_list)
    CALL FillFromTripletList(h_this%data, h_triplet_list%data)
  END SUBROUTINE FillFromTripletList_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Fill in the values of a distributed matrix with the identity matrix.
  PURE SUBROUTINE FillDistributedIdentity_wrp(ih_this) &
       & bind(c,name="FillDistributedIdentity_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    CALL FillDistributedIdentity(h_this%data)
  END SUBROUTINE FillDistributedIdentity_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Fill in the values of a distributed matrix with a permutation matrix.
  SUBROUTINE FillDistributedPermutation_wrp(ih_this, ih_permutation, &
       & permute_rows) bind(c,name="FillDistributedPermutation_wrp")
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: ih_permutation(SIZE_wrp)
    LOGICAL(kind=c_bool), INTENT(IN) :: permute_rows
    TYPE(DistributedSparseMatrix_wrp) :: h_this
    TYPE(Permutation_wrp) :: h_permutation

    h_this = TRANSFER(ih_this,h_this)
    h_permutation = TRANSFER(ih_permutation,h_permutation)

    CALL FillDistributedPermutation(h_this%data,h_permutation%data%index_lookup, &
         & permuterows = LOGICAL(permute_rows))
  END SUBROUTINE FillDistributedPermutation_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the Actual Dimension accessor.
  PURE SUBROUTINE GetActualDimension_wrp(ih_this, mat_dimension) &
       & bind(c,name="GetActualDimension_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(OUT) :: mat_dimension
    TYPE(DistributedSparseMatrix_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    mat_dimension = GetActualDimension(h_this%data)
  END SUBROUTINE GetActualDimension_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Wrap the Logical Dimension accessor.
  PURE SUBROUTINE GetLogicalDimension_wrp(ih_this, mat_dimension) &
       & bind(c,name="GetLogicalDimension_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(OUT) :: mat_dimension
    TYPE(DistributedSparseMatrix_wrp) :: h_this

    h_this = TRANSFER(ih_this,h_this)
    mat_dimension = GetLogicalDimension(h_this%data)
  END SUBROUTINE GetLogicalDimension_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extracts a triplet list of the data that is stored on this process.
  PURE SUBROUTINE GetTripletList_wrp(ih_this, ih_triplet_list) &
       & bind(c,name="GetTripletList_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_triplet_list(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_this
    TYPE(TripletList_wrp) :: h_triplet_list

    h_this = TRANSFER(ih_this,h_this)
    h_triplet_list = TRANSFER(ih_triplet_list,h_triplet_list)
    CALL GetTripletList(h_this%data,h_triplet_list%data)
  END SUBROUTINE GetTripletList_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract an arbitrary block of a matrix into a triplet list.
  SUBROUTINE GetMatrixBlock_wrp(ih_this, ih_triplet_list, start_row, &
       & end_row, start_column, end_column) bind(c,name="GetMatrixBlock_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_this(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_triplet_list(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(IN) :: start_row, end_row
    INTEGER(kind=c_int), INTENT(IN) :: start_column, end_column
    TYPE(DistributedSparseMatrix_wrp) :: h_this
    TYPE(TripletList_wrp) :: h_triplet_list

    h_this = TRANSFER(ih_this,h_this)
    h_triplet_list = TRANSFER(ih_triplet_list,h_triplet_list)
    CALL GetMatrixBlock(h_this%data,h_triplet_list%data, start_row, end_row,&
         & start_column, end_column)
  END SUBROUTINE GetMatrixBlock_wrp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Transpose a sparse matrix.
  SUBROUTINE TransposeDistributedSparseMatrix_wrp(ih_matA,ih_transmat) &
       & bind(c,name="TransposeDistributedSparseMatrix_wrp")
    INTEGER(kind=c_int), INTENT(IN) :: ih_matA(SIZE_wrp)
    INTEGER(kind=c_int), INTENT(INOUT) :: ih_transmat(SIZE_wrp)
    TYPE(DistributedSparseMatrix_wrp) :: h_matA
    TYPE(DistributedSparseMatrix_wrp) :: h_transmat

    h_matA = TRANSFER(ih_matA,h_matA)
    h_transmat = TRANSFER(ih_transmat,h_transmat)
    CALL TransposeDistributedSparseMatrix(h_matA%data,h_transmat%data)
  END SUBROUTINE TransposeDistributedSparseMatrix_wrp
END MODULE DistributedSparseMatrixModule_wrp
