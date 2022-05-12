set -e

if [[ "$TESTOS" == "LINUX" ]]; then
   conda activate ntpoly-conda
fi

cd Build

if [[ "$TESTOS" == "OSX" ]]; then
  cmake -G Ninja .. -DCMAKE_TOOLCHAIN_FILE=../Targets/Mac-python3.cmake \
    -DCMAKE_BUILD_TYPE=Release ;
else
  if  [ ! -z ${NOIALLGATHER+x} ]; then
    cmake -G Ninja .. \
          -DCMAKE_BUILD_TYPE=Debug -DNOIALLGATHER=YES;
  elif [ ! -z ${DEBUG+x} ]; then
    cmake -G Ninja .. -DCMAKE_TOOLCHAIN_FILE=../Targets/Linux.cmake \
          -DCMAKE_BUILD_TYPE=Debug ;
  else
    cmake -G Ninja .. -DCMAKE_TOOLCHAIN_FILE=../Targets/Linux.cmake \
          -DCMAKE_BUILD_TYPE=Release ;
  fi
fi

ninja -v
ninja package
ninja package_source
cd ../
