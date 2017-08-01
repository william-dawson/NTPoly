!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for an iterator that allows one to easily iterate over elements.
MODULE SparseMatrixIteratorModule
  USE DataTypesModule
  USE SparseMatrixModule
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A datatype for iterating over a CSR matrix.
  TYPE, PUBLIC :: SparseMatrixIterator_t
     !! Public Members of the iterator
     !> Total elements to iterate over.
     INTEGER :: total_elements
     !> Current element's row.
     INTEGER :: row
     !> Current element's column.
     INTEGER :: column
     !> Current element's value.
     REAL(NTREAL) :: value
     !! Private Members of the Iterator
     !> Inner counter for iterating.
     INTEGER :: inner_counter
     !> Outer counter for iterating.
     INTEGER :: outer_counter
     !> Counter for iterating over values.
     INTEGER :: total_counter
     !> Elements in the current row.
     INTEGER :: elements_per_inner
  END TYPE SparseMatrixIterator_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: Start
  PUBLIC :: Next
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Associates an iterator with a given matrix and points it to first element.
  !! @param[inout] this the iterator to associate.
  !! @param[in] sparse_matrix the matrix to associate it with.
  PURE SUBROUTINE Start(this,sparse_matrix)
    !! Parameters
    TYPE(SparseMatrixIterator_t), INTENT(inout) :: this
    TYPE(SparseMatrix_t), INTENT(in) :: sparse_matrix

    !! Set Inner Representation
    this%outer_counter = 1
    this%inner_counter = 1
    this%total_counter = 1

    this%elements_per_inner = sparse_matrix%outer_index(this%outer_counter+1) -&
         & sparse_matrix%outer_index(this%outer_counter)

    !! Set Public Values
    this%total_elements = sparse_matrix%outer_index(sparse_matrix%columns+1)
    this%column = this%outer_counter
    this%row = sparse_matrix%inner_index(this%total_counter)
    this%value = sparse_matrix%values(this%total_counter)
  END SUBROUTINE Start
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Gets the next element in a sparse matrix.
  !! @param[inout] this the sparse matrix iterator.
  !! @param[in] sparse_matrix the sparse matrix associated with the iterator.
  !! @todo it's probably a pain to have to pass the sparse matrix here. But
  !! Fortran's requirement that everything that is pointed to has the TARGET
  !! attribute is problematic.
  PURE SUBROUTINE Next(this,sparse_matrix)
    !! Parameters
    TYPE(SparseMatrixIterator_t), INTENT(inout) :: this
    TYPE(SparseMatrix_t), INTENT(in) :: sparse_matrix

    !! Update Inner Representation
    this%total_counter = this%total_counter + 1
    this%inner_counter = this%inner_counter + 1
    IF (this%inner_counter .GT. this%elements_per_inner) THEN
       this%outer_counter = this%outer_counter + 1
       this%elements_per_inner = &
            & sparse_matrix%outer_index(this%outer_counter+1)- &
            & sparse_matrix%outer_index(this%outer_counter)
       this%inner_counter = 1
    END IF

    !! Update Outer Representation
    this%column = this%outer_counter
    this%row = sparse_matrix%inner_index(this%total_counter)
    this%value = sparse_matrix%values(this%total_counter)
  END SUBROUTINE Next
END MODULE SparseMatrixIteratorModule
