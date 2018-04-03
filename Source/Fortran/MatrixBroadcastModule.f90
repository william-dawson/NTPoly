!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Module for broadcasting matrices across processes.
MODULE MatrixBroadcastModule
  USE DataTypesModule
  USE SparseMatrixModule
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: BroadcastMatrix
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Broadcast a matrix from the root process across the communicator.
  !! @param[inout] matrix to broadcast.
  !! @param[in] comm the communicator to broadcast along.
  !! @param[in] root the root process which holds the matrix.
  SUBROUTINE BroadcastMatrix(matrix, comm, root)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(INOUT) :: matrix
    INTEGER, INTENT(INOUT) :: comm
    INTEGER, INTENT(IN) :: root
    !! Local Data
    INTEGER :: rank
    !! Matrix Buffer
    INTEGER :: matrows, matcolumns
    INTEGER :: nnz
    INTEGER :: ierr

    CALL MPI_COMM_RANK(comm, rank, ierr)

    !! Matrix Size
    IF (rank .EQ. root) THEN
       matrows = matrix%rows
       matcolumns = matrix%columns
       nnz = matrix%outer_index(matrix%columns+1)
    END IF
    CALL MPI_Bcast(matrows, 1, MPI_INT, root, comm, ierr)
    CALL MPI_Bcast(matcolumns, 1, MPI_INT, root, comm, ierr)
    CALL MPI_Bcast(nnz, 1, MPI_INT, root, comm, ierr)
    IF (rank .NE. root) THEN
       CALL ConstructZeroSparseMatrix(matrix, matrows, matcolumns)
       DEALLOCATE(matrix%inner_index)
       DEALLOCATE(matrix%values)
       ALLOCATE(matrix%inner_index(nnz))
       ALLOCATE(matrix%values(nnz))
    END IF

    !! Gather Matrix Data
    CALL MPI_Bcast(matrix%outer_index, matcolumns+1, MPI_INT, &
         & root, comm, ierr)
    CALL MPI_Bcast(matrix%inner_index, nnz, MPI_INT, root, comm, ierr)
    CALL MPI_Bcast(matrix%values, nnz, MPINTREAL, root, comm, ierr)
  END SUBROUTINE BroadcastMatrix
END MODULE MatrixBroadcastModule
