if [[ "$TESTOS" == "OSX" ]]; then
  cmake -DCMAKE_TOOLCHAIN_FILE=../Targets/Mac-python3.cmake .. ;
else
  cmake -DCMAKE_TOOLCHAIN_FILE=../Targets/Linux.cmake .. ;
fi
