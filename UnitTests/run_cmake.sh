set -e

if [[ "$TESTOS" == "LINUX" ]]; then
   conda activate ntpoly-conda
fi

cd Build

if [[ "$TESTOS" == "OSX" ]]; then
  cmake -G Ninja .. -DCMAKE_BUILD_TYPE=Release ;
else
  if [[ "${NOIALLGATHER:-0}" -eq 1 ]]; then
    cmake -G Ninja .. -DCMAKE_BUILD_TYPE=Debug -DNOIALLGATHER=YES;
  elif [[ "${DEBUG:-0}" -eq 1 ]]; then
    cmake -G Ninja .. -DCMAKE_TOOLCHAIN_FILE=../Targets/Linux.cmake \
      -DCMAKE_BUILD_TYPE=Debug;
  else
    cmake -G Ninja .. -DCMAKE_TOOLCHAIN_FILE=../Targets/Linux.cmake \
      -DCMAKE_BUILD_TYPE=Release;
  fi
fi

sudo ninja package
cd ../
