  !! Temporary variables
  INTEGER :: alloc_stat
  INTEGER :: num_buckets

  CALL DestructMatrixMemoryPool(this)

  this%columns = columns
  this%rows = rows

  IF (.NOT. PRESENT(sparsity_in)) THEN
     this%hash_size = 1
  ELSE
     this%hash_size = INT(1.0/sparsity_in)
     IF (this%hash_size > columns) this%hash_size = columns
  END IF

  num_buckets = columns/this%hash_size + 1

  !! Allocate
  ALLOCATE(this%pruned_list(columns*rows), stat=alloc_stat)
  ALLOCATE(this%value_array(columns,rows), stat=alloc_stat)
  ALLOCATE(this%dirty_array(columns,rows), stat=alloc_stat)

  ALLOCATE(this%hash_index(columns,rows))
  ALLOCATE(this%inserted_per_bucket(columns,rows))

  this%value_array = 0
  this%hash_index = 0
  this%inserted_per_bucket = 0
  this%dirty_array = .FALSE.
