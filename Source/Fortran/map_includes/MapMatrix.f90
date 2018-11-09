  !! Convert to a triplet list, map the triplet list, fill.
  CALL ConstructEmptyMatrix(outmat, inmat)
  CALL GetMatrixTripletList(inmat, inlist)
#ifdef MAPARRAY
  CALL MapTripletList(inlist, outlist, proc, supp_in=supp_in, &
       & num_slices_in=inmat%process_grid%num_process_slices, &
       & my_slice_in=inmat%process_grid%my_slice)
#else
  CALL MapTripletList(inlist, outlist, proc, &
       & num_slices_in=inmat%process_grid%num_process_slices, &
       & my_slice_in=inmat%process_grid%my_slice)
#endif
  CALL FillMatrixFromTripletList(outmat, outlist, preduplicated_in=.FALSE.)

  !! Cleanup
  CALL DestructTripletList(inlist)
  CALL DestructTripletList(outlist)
