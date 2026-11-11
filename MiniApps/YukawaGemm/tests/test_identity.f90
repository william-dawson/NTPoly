PROGRAM IdentityRegression
  USE NTMPIModule
  USE DataTypesModule, ONLY: NTREAL
  USE ProcessGridModule
  USE PSMatrixModule
  USE TripletListModule
  IMPLICIT NONE
  TYPE(Matrix_ps) :: a
  TYPE(TripletList_r) :: real_entries
  TYPE(TripletList_c) :: complex_entries
  INTEGER :: ierr,nr,rank,orientation,pr,pc,n,kind,j,local_count,total
  INTEGER, ALLOCATABLE :: diagonal_counts(:)
  INTEGER, PARAMETER :: sizes(3)=[1,5,8]
  CALL MPI_Init(ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD,nr,ierr)
  DO orientation=1,2
    pr=1; pc=nr
    IF(orientation==2) THEN
      pr=nr; pc=1
    END IF
    CALL ConstructProcessGrid(MPI_COMM_WORLD,pr,pc,1)
    DO n=1,SIZE(sizes)
      DO kind=1,2
        ALLOCATE(diagonal_counts(sizes(n))); diagonal_counts=0
        CALL ConstructEmptyMatrix(a,sizes(n),is_complex_in=kind==2)
        CALL FillMatrixIdentity(a)
        IF(kind==1) THEN
          CALL GetMatrixTripletList(a,real_entries)
          local_count=real_entries%CurrentSize
          DO j=1,local_count
            IF(real_entries%DATA(j)%index_row/=real_entries%DATA(j)%index_column .OR. &
               real_entries%DATA(j)%index_row<1 .OR. real_entries%DATA(j)%index_row>sizes(n) .OR. &
               real_entries%DATA(j)%point_value/=1.0_NTREAL) &
               CALL MPI_Abort(MPI_COMM_WORLD,1,ierr)
            diagonal_counts(real_entries%DATA(j)%index_row)= &
              diagonal_counts(real_entries%DATA(j)%index_row)+1
          END DO
          CALL DestructTripletList(real_entries)
        ELSE
          CALL GetMatrixTripletList(a,complex_entries)
          local_count=complex_entries%CurrentSize
          DO j=1,local_count
            IF(complex_entries%DATA(j)%index_row/=complex_entries%DATA(j)%index_column .OR. &
               complex_entries%DATA(j)%index_row<1 .OR. complex_entries%DATA(j)%index_row>sizes(n) .OR. &
               complex_entries%DATA(j)%point_value/=1.0_NTREAL) &
               CALL MPI_Abort(MPI_COMM_WORLD,2,ierr)
            diagonal_counts(complex_entries%DATA(j)%index_row)= &
              diagonal_counts(complex_entries%DATA(j)%index_row)+1
          END DO
          CALL DestructTripletList(complex_entries)
        END IF
        CALL MPI_Allreduce(local_count,total,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        IF(total/=sizes(n)) CALL MPI_Abort(MPI_COMM_WORLD,3,ierr)
        CALL MPI_Allreduce(MPI_IN_PLACE,diagonal_counts,sizes(n),MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        IF(ANY(diagonal_counts/=1)) CALL MPI_Abort(MPI_COMM_WORLD,4,ierr)
        DEALLOCATE(diagonal_counts)
        CALL DestructMatrix(a)
      END DO
    END DO
    CALL DestructProcessGrid()
  END DO
  IF(rank==0) PRINT *, 'Identity real/complex, padded and empty local blocks: PASS',nr
  CALL MPI_Finalize(ierr)
END PROGRAM
