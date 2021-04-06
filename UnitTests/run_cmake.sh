if [[ "$TESTOS" == "OSX" ]]; then
  cmake -G ninja .. -DCMAKE_TOOLCHAIN_FILE=../Targets/Mac-python3.cmake \
    -DCMAKE_BUILD_TYPE=Release ;
else
  if  [ ! -z ${NOIALLGATHER+x} ]; then
    cmake -G ninja .. \
          -DCMAKE_BUILD_TYPE=Debug -DNOIALLGATHER=YES;
  elif [ ! -z ${DEBUG+x} ]; then
    cmake -G ninja .. -DCMAKE_TOOLCHAIN_FILE=../Targets/Linux.cmake \
          -DCMAKE_BUILD_TYPE=Debug ;
  else
    cmake -G ninja .. -DCMAKE_TOOLCHAIN_FILE=../Targets/Linux.cmake \
          -DCMAKE_BUILD_TYPE=Release ;
  fi
fi
