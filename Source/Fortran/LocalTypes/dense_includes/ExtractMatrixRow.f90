  CALL row_out%InitEmpty(1, this%columns)
  row_out%data(1,:) = this%data(row_number, :)
