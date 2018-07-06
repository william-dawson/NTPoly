  !! Helper variables
  INTEGER :: outer_counter, inner_counter
  INTEGER :: elements_per_inner
  INTEGER :: total_counter
  INTEGER :: size_of_this

  size_of_this = this%outer_index(this%columns+1)
  CALL ConstructTripletList(triplet_list, size_of_this)

  total_counter = 1
  DO outer_counter = 1, this%columns
     elements_per_inner = this%outer_index(outer_counter+1) - &
          & this%outer_index(outer_counter)
     DO inner_counter = 1, elements_per_inner
        temporary%index_column = outer_counter
        temporary%index_row = this%inner_index(total_counter)
        temporary%point_value = this%values(total_counter)
        triplet_list%data(total_counter) = temporary
        total_counter = total_counter + 1
     END DO
  END DO
