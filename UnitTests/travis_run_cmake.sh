if [[ "$TESTOS" == "OSX" ]]; then
  cmake -DCMAKE_TOOLCHAIN_FILE=../Targets/Mac-python2.cmake .. ;
else
  if  [ ! -z ${NOIALLGATHER+x} ]; then
    cmake -DCMAKE_TOOLCHAIN_FILE=../Targets/Linux.cmake \
          -DNOIALLGATHER=YES .. ;
  elif [ -z ${FORTRAN_ONLY+x} ]; then
    cmake -DCMAKE_TOOLCHAIN_FILE=../Targets/Linux.cmake .. ;
  else
    cmake -DCMAKE_TOOLCHAIN_FILE=../Targets/Linux.cmake \
          -DFORTRAN_ONLY=YES .. ;
  fi
fi
