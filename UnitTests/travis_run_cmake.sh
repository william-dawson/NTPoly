if [[ "$TESTOS" == "OSX" ]]; then
  cmake -DCMAKE_TOOLCHAIN_FILE=../Targets/Mac-python2.cmake .. ;
else
  cmake -DCMAKE_TOOLCHAIN_FILE=../Targets/Linux.cmake .. ;
fi
