CALL column_out%InitEmpty(this%rows, 1)
column_out%data(:,1) = this%data(:, column_number)
