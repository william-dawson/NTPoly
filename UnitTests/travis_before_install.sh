
if [[ "$TESTOS" == "LINUX" ]]; then
  sudo apt-get install gfortran
  sudo apt-get install libblas-dev liblapack-dev
  if [[ "$MPICH" == "1" ]]; then
    sudo apt-get install mpich libhwloc-plugins libmpich-dev
  else
    sudo apt-get install openmpi-bin libopenmpi-dev libhwloc-contrib-plugins openmpi-doc opencl-icd
  fi
fi

if [[ "$TESTOS" == "OSX" ]]; then
  brew install gcc
  brew link --overwrite gcc
  brew install open-mpi
  brew install doxygen
  brew install cmake
  brew install swig
  sudo pip2 install scipy --upgrade --no-cache-dir
  sudo pip2 install numpy --upgrade --no-cache-dir
  sudo pip2 install mpi4py --upgrade --no-cache-dir
  sudo pip2 install flake8 --upgrade --no-cache-dir
else
  sudo ldconfig
  sudo apt-get install python-dev python-pip python-all-dev python-setuptools python-wheel
  sudo apt-get install swig
  sudo pip install scipy --upgrade
  sudo pip install numpy --upgrade
  sudo pip install mpi4py --upgrade
  sudo pip install flake8 --upgrade
fi

test -n $CC  && unset CC
test -n $FCC  && unset FCC
test -n $CXX && unset CXX
