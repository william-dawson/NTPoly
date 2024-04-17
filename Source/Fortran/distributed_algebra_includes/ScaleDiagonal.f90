INTEGER :: II, row

!! Merge to the local block
CALL MergeMatrixLocalBlocks(this, lmat)

!! Filter out the triplets that aren't stored locally
CALL ConstructTripletList(filtered)
DO II = 1, tlist%CurrentSize
    CALL GetTripletAt(tlist, II, trip)
    row = trip%index_row
    IF (row .GE. this%start_row .AND. row .LT. this%end_row) THEN
        trip%index_row = trip%index_row - this%start_row + 1
        trip%index_column = trip%index_row
        CALL AppendToTripletList(filtered, trip)
    END IF
END DO

!! Scale
CALL MatrixDiagonalScale(lmat, filtered)

!! Split
CALL SplitMatrixToLocalBlocks(this, lmat)
CALL DestructMatrix(lmat)
CALL DestructTripletList(filtered)