  tripB%CurrentSize = tripA%CurrentSize

  !! We only will allocate as much space as needed, and not the additional
  !! buffer.
  ALLOCATE(tripB%DATA(tripB%CurrentSize))
  tripB%DATA(:tripB%CurrentSize) = tripA%DATA(:tripB%CurrentSize)