  !! Local Variables
  INTEGER :: num_buckets

  this%hash_size = INT(1.0 / sparsity)
  IF (this%hash_size .GT. this%columns) this%hash_size = this%columns
  num_buckets = this%columns / this%hash_size + 1
