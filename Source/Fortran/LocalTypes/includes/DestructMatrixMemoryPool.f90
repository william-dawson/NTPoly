  !! Perform deallocations.
  IF (ALLOCATED(this%pruned_list)) DEALLOCATE(this%pruned_list)
  IF (ALLOCATED(this%value_array)) DEALLOCATE(this%value_array)
  IF (ALLOCATED(this%dirty_array)) DEALLOCATE(this%dirty_array)
  IF (ALLOCATED(this%hash_index)) DEALLOCATE(this%hash_index)
  IF (ALLOCATED(this%inserted_per_bucket)) DEALLOCATE(this%inserted_per_bucket)
