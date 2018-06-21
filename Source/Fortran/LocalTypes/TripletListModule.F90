!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Storing Lists of Triplets.
MODULE TripletListModule
  USE DataTypesModule, ONLY: NTREAL, MPINTREAL, NTCOMPLEX, MPINTCOMPLEX
  USE TripletModule, ONLY : Triplet, Triplet_r, Triplet_c
  USE MatrixMarketModule, ONLY : MM_SYMMETRIC, MM_SKEW_SYMMETRIC, MM_HERMITIAN
  USE TimerModule, ONLY : StartTimer, StopTimer
  USE ISO_C_BINDING, ONLY : c_int
  USE MPI
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE, PUBLIC :: TripletList
     !> Internal representation of the data (real).
     TYPE(Triplet_r), DIMENSION(:), ALLOCATABLE :: data_r
     !> Internal representation of the data (complex).
     TYPE(Triplet_c), DIMENSION(:), ALLOCATABLE :: data_c
     !> Current number of elements in the triplet list
     INTEGER :: CurrentSize
     !> Whether we're storing real or complex data.
     LOGICAL :: IsReal
   CONTAINS
     PROCEDURE :: Init => ConstructTripletList
     PROCEDURE :: Destruct => DestructTripletList
     PROCEDURE :: Resize => ResizeTripletList
     PROCEDURE :: Append => AppendToTripletList
     PROCEDURE :: Set => SetTripletAt
     PROCEDURE :: Get => GetTripletAt
     PROCEDURE, NOPASS :: Sort => SortTripletList
     PROCEDURE :: Symmetrize => SymmetrizeTripletList
     PROCEDURE :: GetSize => GetTripletListSize
     PROCEDURE, NOPASS :: Redistribute => RedistributeTripletLists
     PROCEDURE :: Shift => ShiftTripletList
  END TYPE TripletList
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a triplet list.
  !! @param[inout] this the triplet list to construct.
  !! @param[in] size_in the length of the triplet list (optional, default=0).
  !! @param[in] is_real whether or not the triplet list stores real data.
  PURE SUBROUTINE ConstructTripletList(this, size_in, is_real)
    !! Parameters
    CLASS(TripletList), INTENT(INOUT) :: this
    INTEGER(kind=c_int), INTENT(IN), OPTIONAL :: size_in
    LOGICAL, INTENT(IN), OPTIONAL :: is_real

    CALL this%Destruct

    IF (PRESENT(size_in)) THEN
       this%CurrentSize  = size_in
    ELSE
       this%CurrentSize  = 0
    END IF
    IF (PRESENT(is_real)) THEN
       this%IsReal = is_real
    ELSE
       this%IsReal = .TRUE.
    END IF

    IF (this%IsReal) THEN
       ALLOCATE(this%data_r(this%CurrentSize))
    ELSE
       ALLOCATE(this%data_c(this%CurrentSize))
    END IF

  END SUBROUTINE ConstructTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destructs a triplet list.
  !! @param[inout] this the triplet list to destruct.
  PURE SUBROUTINE DestructTripletList(this)
    !! Parameters
    CLASS(TripletList), INTENT(INOUT) :: this

    IF (ALLOCATED(this%data_r)) DEALLOCATE(this%data_r)
    IF (ALLOCATED(this%data_c)) DEALLOCATE(this%data_c)
    this%CurrentSize = 0

  END SUBROUTINE DestructTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Increase the size of a triplet list.
  !! @param[inout] this the triplet list to resize.
  !! @param[in] size to resize to.
  PURE SUBROUTINE ResizeTripletList(this, size)
    !! Parameters
    CLASS(TripletList), INTENT(INOUT) :: this
    INTEGER(KIND=c_int), INTENT(IN) :: size
    !! Local Data
    TYPE(Triplet_r), DIMENSION(:), ALLOCATABLE :: temporary_data_r
    TYPE(Triplet_c), DIMENSION(:), ALLOCATABLE :: temporary_data_c

    IF (this%IsReal) THEN
#define temporary_data temporary_data_r
#define DATA data_r
#include "triplet_includes/ResizeTripletList.f90"
#undef DATA
#undef temporary_data
    ELSE
#define temporary_data temporary_data_c
#define DATA data_c
#include "triplet_includes/ResizeTripletList.f90"
#undef DATA
#undef temporary_data
    END IF

  END SUBROUTINE ResizeTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Add a value to the end of the triplet list.
  !! @param[inout] this the triplet list to append to.
  !! @param[in] triplet_value the value to append.
  PURE SUBROUTINE AppendToTripletList(this, triplet_value)
    !! Parameters
    CLASS(TripletList), INTENT(INOUT) :: this
    CLASS(Triplet), INTENT(IN)         :: triplet_value
    !! Local data
    INTEGER :: new_size
    INTEGER :: old_size

    IF (this%IsReal) THEN
       old_size = SIZE(this%data_r)
    ELSE
       old_size = SIZE(this%data_c)
    END IF

    !! First, check if we need to allocate more memory
    IF (this%CurrentSize+1 .GT. old_size) THEN
       IF (old_size .EQ. 0) THEN
          new_size = 1
       ELSE IF (old_size .EQ. 1) THEN
          new_size = 2
       ELSE
          new_size = INT(old_size*1.5)
       END IF
       CALL this%Resize(new_size)
    END IF

    !! Append
    this%CurrentSize = this%CurrentSize+1

    SELECT TYPE(triplet_value)
    CLASS IS (Triplet_r)
      this%data_r(this%CurrentSize) = triplet_value
    CLASS IS (Triplet_c)
      this%data_c(this%CurrentSize) = triplet_value
    END SELECT

  END SUBROUTINE AppendToTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of a triplet at a particular index.
  !! @param[inout] this the triplet list to set.
  !! @param[in] index the index at which to set the triplet.
  !! @param[in] triplet_value the value of the triplet to set.
  PURE SUBROUTINE SetTripletAt(this,index,triplet_value)
    !! Parameters
    CLASS(TripletList), INTENT(INOUT) :: this
    INTEGER(KIND=c_int), INTENT(IN)    :: index
    CLASS(Triplet), INTENT(IN)        :: triplet_value

    SELECT TYPE(triplet_value)
    CLASS IS (Triplet_r)
      this%data_r(index) = triplet_value
    CLASS IS (Triplet_c)
      this%data_c(index) = triplet_value
    END SELECT

  END SUBROUTINE SetTripletAt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the value of a triplet at a particular index.
  !! @param[in] this the triplet list to get the value from.
  !! @param[in] index the index from which to get the triplet.
  !! @param[out] triplet_value the extracted triplet value.
  PURE SUBROUTINE GetTripletAt(this,index,triplet_value)
    !! Parameters
    CLASS(TripletList), INTENT(IN) :: this
    INTEGER(kind=c_int), INTENT(IN) :: index
    CLASS(Triplet), INTENT(INOUT)   :: triplet_value

    SELECT TYPE(triplet_value)
    CLASS IS (Triplet_r)
      triplet_value%index_row = this%data_r(index)%index_row
      triplet_value%index_column = this%data_r(index)%index_column
      triplet_value%point_value = this%data_r(index)%point_value
    CLASS IS (Triplet_c)
      triplet_value%index_row = this%data_c(index)%index_row
      triplet_value%index_column = this%data_c(index)%index_column
      triplet_value%point_value = this%data_c(index)%point_value
    END SELECT

  END SUBROUTINE GetTripletAt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sorts a triplet list by index values.
  !! Implementation is based on bucket sort. This is why it needs the number of
  !! matrix columns. Bubble sort is used within a bucket.
  !! @param[in] input_list list to be sorted.
  !! @param[in] matrix_columns this is the highest column value in the list.
  !! @param[in] matrix_rows this is the highest row value in the list.
  !! @param[in] bubble_in false if you don't need the final bubble sort.
  !! @param[out] sorted_list a now sorted version of the list. This routine
  !! will allocate it.
  PURE SUBROUTINE SortTripletList(input_list, matrix_columns, matrix_rows, &
       & sorted_list, bubble_in)
    !! Parameters
    CLASS(TripletList), INTENT(IN)  :: input_list
    INTEGER, INTENT(IN) :: matrix_columns
    INTEGER, INTENT(IN) :: matrix_rows
    TYPE(TripletList), INTENT(OUT) :: sorted_list
    LOGICAL, OPTIONAL, INTENT(IN) :: bubble_in
    !! Local Data
    TYPE(Triplet_r) :: temporary_r
    TYPE(Triplet_c) :: temporary_c
    !! Local Data
    LOGICAL :: bubble
    LOGICAL :: swap_occured
    INTEGER, DIMENSION(:), ALLOCATABLE :: values_per_row
    INTEGER, DIMENSION(:), ALLOCATABLE :: offset_array
    INTEGER, DIMENSION(:), ALLOCATABLE :: inserted_per_row
    !! Counters and temporary variables
    INTEGER :: counter
    INTEGER :: temp_index
    INTEGER :: alloc_stat
    INTEGER :: list_length

    IF (PRESENT(bubble_in)) THEN
       bubble = bubble_in
    ELSE
       bubble = .TRUE.
    END IF

    list_length = input_list%CurrentSize

    IF (input_list%IsReal) THEN
#define DATA data_r
#define temporary temporary_r
#include "triplet_includes/SortTripletList.f90"
#undef temporary
#undef DATA
    ELSE
#define DATA data_c
#define temporary temporary_c
#include "triplet_includes/SortTripletList.f90"
#undef temporary
#undef DATA
    END IF

  END SUBROUTINE SortTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sort a triplet list assuming that the matrix it corresponds to is nearly
  !! dense.
  !! @param[in] input_list the list
  !! @param[in] matrix_columns for the corresponding matrix.
  !! @param[in] matrix_rows for the corresponding matrix.
  !! @param[out] sorted_list sorted and ready to use for building matrices.
  PURE SUBROUTINE SortDenseTripletList(input_list, matrix_columns, &
       & matrix_rows, sorted_list)
    !! Parameters
    TYPE(TripletList), INTENT(IN)  :: input_list
    INTEGER, INTENT(IN) :: matrix_columns
    INTEGER, INTENT(IN) :: matrix_rows
    TYPE(TripletList), INTENT(OUT) :: sorted_list
    !! Local Variables
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: value_buffer_r
    COMPLEX(NTCOMPLEX), DIMENSION(:,:), ALLOCATABLE :: value_buffer_c
    !! Local Data
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: dirty_buffer
    INTEGER :: list_length
    INTEGER :: row, col, ind
    INTEGER :: II, JJ

    IF (input_list%IsReal) THEN
#define DATA data_r
#define value_buffer value_buffer_r
#include "triplet_includes/SortDenseTripletList.f90"
#undef value_buffer
#undef DATA
    ELSE
#define DATA data_c
#define value_buffer value_buffer_c
#include "triplet_includes/SortDenseTripletList.f90"
#undef value_buffer
#undef DATA
    END IF

  END SUBROUTINE SortDenseTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Symmetrizes an unsymmetric triplet list according to the specified
  !! symmetry type.
  !! @param[inout] triplet_list list to be symmetrized.
  !! @param[in] pattern_type type of symmetry.
  SUBROUTINE SymmetrizeTripletList(this, pattern_type)
    !! Parameters
    CLASS(TripletList), INTENT(INOUT)  :: this
    INTEGER, INTENT(IN) :: pattern_type
    !! Local variables
    CLASS(Triplet), POINTER :: temporary
    TYPE(Triplet_r), TARGET :: temporary_r
    TYPE(Triplet_c), TARGET :: temporary_c
    INTEGER :: counter
    INTEGER :: initial_size

    initial_size = this%CurrentSize
    DO counter = 1, initial_size
       IF (this%IsReal) THEN
          CALL this%Get(counter,temporary_r)
          temporary => temporary_r
       ELSE
          CALL this%Get(counter,temporary_c)
          temporary => temporary_c
       END IF
       IF (temporary%index_row .EQ. temporary%index_column) CONTINUE
       CALL temporary%Transpose()

       SELECT CASE(pattern_type)
       CASE(MM_HERMITIAN)
           CALL temporary%Conjg()
       CASE(MM_SKEW_SYMMETRIC)
           CALL temporary%Scale(REAL(-1.0,NTREAL))
       END SELECT

       CALL this%Append(temporary)
    END DO
  END SUBROUTINE SymmetrizeTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of entries in a triplet list.
  !! @param[in] triplet_list list to get the size of.
  !! @return list_size the number of entries in the triplet list.
  PURE FUNCTION GetTripletListSize(triplet_list) RESULT(list_size)
    !! Parameters
    CLASS(TripletList), INTENT(IN)  :: triplet_list
    INTEGER :: list_size

    list_size = triplet_list%CurrentSize

  END FUNCTION GetTripletListSize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Redistribute some triplet lists amongst a set of processors.
  !! Takes in a list of triplet lists, one list for each processor. Then the
  !! all to all redistribution is performed along the given communicator.
  !! @param[in] triplet_lists a list of triplet lists, one for each process.
  !! @param[inout] comm the mpi communicator to redistribute along.
  !! @param[out] local_data_out the resulting local triplet list.
  SUBROUTINE RedistributeTripletLists(triplet_lists, comm, local_data_out)
    !! Parameters
    CLASS(TripletList), DIMENSION(:), INTENT(IN) :: triplet_lists
    INTEGER, INTENT(INOUT) :: comm
    CLASS(TripletList), INTENT(INOUT) :: local_data_out
    !! Local data (type specific)
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: send_buffer_val
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: recv_buffer_val
    TYPE(Triplet_r) :: temp_triplet_r
    TYPE(Triplet_r) :: temp_triplet_c
    !! Local Data - Offsets
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_per_process
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_offsets
    INTEGER, DIMENSION(:), ALLOCATABLE :: recv_per_process
    INTEGER, DIMENSION(:), ALLOCATABLE :: recv_offsets
    !! Local Data - Send/Recv Buffers
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_buffer_row
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_buffer_col
    INTEGER, DIMENSION(:), ALLOCATABLE :: recv_buffer_row
    INTEGER, DIMENSION(:), ALLOCATABLE :: recv_buffer_col
    !! ETC
    INTEGER :: num_processes
    INTEGER :: counter, inner_counter, insert_pt
    INTEGER :: mpi_error
    LOGICAL :: is_real

    IF (SIZE(triplet_lists) .GT. 0) THEN
      is_real = triplet_lists(1)%IsReal
    ELSE
      is_real = .TRUE.
    END IF

    !! Allocate Size Buffers
    CALL MPI_COMM_SIZE(comm, num_processes, mpi_error)
    ALLOCATE(send_per_process(num_processes))
    ALLOCATE(send_offsets(num_processes))
    ALLOCATE(recv_per_process(num_processes))
    ALLOCATE(recv_offsets(num_processes))

    !! Figure Out How Much Data Gets Sent
    DO counter = 1, num_processes
       send_per_process(counter) = triplet_lists(counter)%CurrentSize
    END DO
    send_offsets(1) = 0
    DO counter = 2, num_processes
       send_offsets(counter) = send_offsets(counter-1) + &
            & send_per_process(counter-1)
    END DO

    !! Figure Out How Much Data Gets Received
    CALL MPI_ALLTOALL(send_per_process, 1, MPI_INT, recv_per_process, 1, &
         & MPI_INT, comm, mpi_error)
    recv_offsets(1) = 0
    DO counter = 2, num_processes
       recv_offsets(counter) = recv_offsets(counter-1) + &
            & recv_per_process(counter-1)
    END DO

    !! Allocate And Fill Send Buffers
    ALLOCATE(send_buffer_row(SUM(send_per_process)))
    ALLOCATE(send_buffer_col(SUM(send_per_process)))
    ALLOCATE(send_buffer_val(SUM(send_per_process)))
    ALLOCATE(recv_buffer_row(SUM(recv_per_process)))
    ALLOCATE(recv_buffer_col(SUM(recv_per_process)))
    ALLOCATE(recv_buffer_val(SUM(recv_per_process)))

    IF (is_real) THEN
#define DATA data_r
#define temp_triplet temp_triplet_r
#define MPIDATATYPE MPINTREAL
#include "triplet_includes/RedistributeTripletLists.f90"
#undef temp_triplet
#undef MPIDATATYPE
#undef DATA
    ELSE
#define DATA data_c
#define temp_triplet temp_triplet_c
#define MPIDATATYPE MPINTCOMPLEX
#include "triplet_includes/RedistributeTripletLists.f90"
#undef temp_triplet
#undef MPIDATATYPE
#undef DATA
    END IF

    !! Cleanup
    DEALLOCATE(send_per_process)
    DEALLOCATE(send_offsets)
    DEALLOCATE(recv_per_process)
    DEALLOCATE(recv_offsets)
    DEALLOCATE(send_buffer_row)
    DEALLOCATE(send_buffer_col)
    DEALLOCATE(send_buffer_val)
    DEALLOCATE(recv_buffer_row)
    DEALLOCATE(recv_buffer_col)
    DEALLOCATE(recv_buffer_val)

  END SUBROUTINE RedistributeTripletLists
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Shift the rows and columns of a triplet list by set values.
  !! Frequently, we have a triplet list that comes from the global matrix which
  !! we would like to shift into a local matrix. In that case, just pass
  !! the negative of the starting row and column (plus 1) to this routine.
  PURE SUBROUTINE ShiftTripletList(this, row_shift, column_shift)
    !! Parameters
    CLASS(TripletList), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: row_shift
    INTEGER, INTENT(IN) :: column_shift
    !! Local Variables
    INTEGER :: II

    !! Loop
    IF (this%IsReal) THEN
      DO II = 1, this%CurrentSize
         this%data_r(II)%index_row = this%data_r(II)%index_row + row_shift
         this%data_r(II)%index_column = &
              this%data_r(II)%index_column + column_shift
      END DO
    ELSE
      DO II = 1, this%CurrentSize
         this%data_c(II)%index_row = this%data_c(II)%index_row + row_shift
         this%data_c(II)%index_column = &
              this%data_c(II)%index_column + column_shift
      END DO
    END IF

  END SUBROUTINE ShiftTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE TripletListModule
