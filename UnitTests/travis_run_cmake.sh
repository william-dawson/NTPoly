if [[ "$TESTOS" == "OSX" ]]; then
  cmake -DCMAKE_TOOLCHAIN_FILE=../Targets/Mac-python2.cmake .. ;
else
  if [ -z ${FORTRAN_ONLY+x} ]; then
    cmake -DCMAKE_TOOLCHAIN_FILE=../Targets/Linux.cmake \
          -DFORTRAN_ONLY=YES.. ;
  else
    cmake -DCMAKE_TOOLCHAIN_FILE=../Targets/Linux.cmake .. ;
  fi
fi
