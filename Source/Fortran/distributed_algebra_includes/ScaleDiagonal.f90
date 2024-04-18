INTEGER :: II, col

!! Merge to the local block
CALL MergeMatrixLocalBlocks(this, lmat)

!! Filter out the triplets that aren't stored locally
CALL ConstructTripletList(filtered)
DO II = 1, tlist%CurrentSize
    CALL GetTripletAt(tlist, II, trip)
    col = trip%index_column
    IF (col .GE. this%start_column .AND. col .LT. this%end_column) THEN
        trip%index_column = trip%index_column - this%start_column + 1
        trip%index_row = trip%index_column
        CALL AppendToTripletList(filtered, trip)
    END IF
END DO

!! Scale
CALL MatrixDiagonalScale(lmat, filtered)

!! Split
CALL SplitMatrixToLocalBlocks(this, lmat)
CALL DestructMatrix(lmat)
CALL DestructTripletList(filtered)
