PROGRAM ResizeGridTest
  USE NTMPIModule
  USE ProcessGridModule
  USE SMatrixModule, ONLY: Matrix_lsr,Matrix_lsc,DestructLocal=>DestructMatrix
  USE PSMatrixModule
  USE DataTypesModule, ONLY: NTREAL
  USE TripletListModule
  USE TripletModule
  IMPLICIT NONE
  TYPE(ProcessGrid_t), TARGET :: grid
  TYPE(Matrix_ps) :: matrix
  TYPE(Matrix_lsr) :: local_r
  TYPE(Matrix_lsc) :: local_c
  TYPE(TripletList_c) :: list
  TYPE(Triplet_c) :: entry
  INTEGER :: ierr,nr,rank,kind,i,target,stage
  CALL MPI_Init(ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD,nr,ierr);CALL MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  ! Deliberately leave the default grid unconstructed.
  CALL ConstructNewProcessGrid(grid,MPI_COMM_WORLD,1,nr,1)
  DO kind=0,1
    CALL ConstructEmptyMatrix(matrix,8,process_grid_in=grid,is_complex_in=kind==1)
    CALL FillMatrixIdentity(matrix)
    DO stage=1,2
      target=5
      IF(stage==2) target=11
      CALL ResizeMatrix(matrix,target)
      IF(.NOT.ASSOCIATED(matrix%process_grid,grid)) CALL MPI_Abort(MPI_COMM_WORLD,1,ierr)
      IF(matrix%is_complex .NEQV. (kind==1)) CALL MPI_Abort(MPI_COMM_WORLD,2,ierr)
      IF(matrix%actual_matrix_dimension/=target) CALL MPI_Abort(MPI_COMM_WORLD,3,ierr)
      IF(GetMatrixSize(matrix)/=5) CALL MPI_Abort(MPI_COMM_WORLD,4,ierr)
      ! Replacement preserves values, and the destination owns its arrays.
      IF(kind==0) THEN
        CALL MergeMatrixLocalBlocks(matrix,local_r)
        CALL SplitMatrixToLocalBlocks(matrix,local_r)
        local_r%values=-7.0_NTREAL
        CALL DestructLocal(local_r)
      ELSE
        CALL MergeMatrixLocalBlocks(matrix,local_c)
        CALL SplitMatrixToLocalBlocks(matrix,local_c)
        local_c%values=CMPLX(-7.0_NTREAL,.25_NTREAL,NTREAL)
        CALL DestructLocal(local_c)
      END IF
      CALL GetMatrixTripletList(matrix,list)
      DO i=1,list%CurrentSize
        entry=list%DATA(i)
        IF(entry%index_row/=entry%index_column .OR. entry%index_row>5 .OR. &
           ABS(entry%point_value-1.0_NTREAL)>1e-15_NTREAL) CALL MPI_Abort(MPI_COMM_WORLD,5,ierr)
      END DO
      CALL DestructTripletList(list)
    END DO
    CALL DestructMatrix(matrix)
  END DO
  CALL DestructProcessGrid(grid)
  IF(rank==0) PRINT *, 'RESIZE CUSTOM GRID REAL/COMPLEX SHRINK/GROW PASS; SPLIT OWNERSHIP PASS'
  CALL MPI_Finalize(ierr)
END PROGRAM
