  !! Local Data
  INTEGER :: II, JJ
  INTEGER :: elements_per_inner

  !! Allocate Space For Result
  norm_per_column = 0

  !! Iterate Over Local Data
  DO II = 1, this%columns
     elements_per_inner = this%outer_index(II + 1) - this%outer_index(II)
     DO JJ = 1, elements_per_inner
        temp_value = this%values(this%outer_index(II) + JJ)
        norm_per_column(II) = norm_per_column(II) + ABS(temp_value)
     END DO
  END DO
