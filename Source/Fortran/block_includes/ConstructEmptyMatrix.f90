    CALL DestructMatrix(this)

    IF (base_case) THEN
       ALLOCATE(this%base_data(block_rows, block_cols))
    ELSE
       ALLOCATE(this%h_data(block_rows, block_cols))
    END IF