  !! Helper variables
  INTEGER :: II, JJ, KK
  INTEGER :: elements_per_inner
  INTEGER :: size_of_this

  size_of_this = this%outer_index(this%columns + 1)
  CALL ConstructTripletList(triplet_list, size_of_this)

  KK = 1
  DO II = 1, this%columns
     elements_per_inner = this%outer_index(II + 1) - this%outer_index(II)
     DO JJ = 1, elements_per_inner
        temporary%index_column = II
        temporary%index_row = this%inner_index(KK)
        temporary%point_value = this%values(KK)
        triplet_list%DATA(KK) = temporary
        KK = KK + 1
     END DO
  END DO
