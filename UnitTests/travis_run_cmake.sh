if [[ "$TESTOS" == "OSX" ]]; then
  cmake .. -DCMAKE_TOOLCHAIN_FILE=../Targets/Mac-python2.cmake \
    -DCMAKE_BUILD_TYPE=Release ;
else
  if  [ ! -z ${NOIALLGATHER+x} ]; then
    cmake .. -DCMAKE_TOOLCHAIN_FILE=../Targets/Linux.cmake \
          -DNOIALLGATHER=YES -DCMAKE_BUILD_TYPE=Debug ;
  elif [ -z ${FORTRAN_ONLY+x} ]; then
    cmake .. -DCMAKE_TOOLCHAIN_FILE=../Targets/Linux.cmake \
          -DCMAKE_BUILD_TYPE=Release ;
  else
    cmake .. -DCMAKE_TOOLCHAIN_FILE=../Targets/Linux.cmake \
          -DFORTRAN_ONLY=YES -DCMAKE_BUILD_TYPE=Release ;
  fi
fi
