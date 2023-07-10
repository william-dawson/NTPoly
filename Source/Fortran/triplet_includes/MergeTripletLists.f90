  INTEGER :: new_size
  INTEGER :: II, JJ, KK

  !! Figure out how big the list should be
  new_size = 0
  DO II = 1, SIZE(tlists)
     new_size = new_size + tlists(II)%CurrentSize
  END DO

  !! Allocate
  CALL ConstructTripletList(this, new_size)

  !! Copy
  KK = 1
  JJ = 1
  DO II = 1, SIZE(tlists)
     DO JJ = 1, tlists(II)%CurrentSize
        this%DATA(KK) = tlists(II)%DATA(JJ)
        KK = KK + 1
     END DO
  END DO
