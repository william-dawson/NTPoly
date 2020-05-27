
if [[ "$TESTOS" == "LINUX" ]]; then
  sudo apt-get update
  sudo apt-get install gfortran
  sudo apt-get install libblas-dev liblapack-dev
  sudo apt-get install gawk
  sudo apt-get install clang-format
  sudo apt-get install emacs
  if [[ "$MPICH" == "1" ]]; then
    sudo apt-get install mpich libmpich-dev
  else
    sudo apt-get install openmpi-bin libopenmpi-dev
  fi
fi

if [[ "$TESTOS" == "OSX" ]]; then
  brew install gcc
  brew link --overwrite gcc
  brew install open-mpi
  brew install doxygen
  brew install cmake
  brew install swig
  brew install clang-format
  brew install emacs
  sudo pip2 install numpy --upgrade --no-cache-dir
  sudo pip2 install scipy --upgrade --no-cache-dir
  sudo pip2 install mpi4py --upgrade --no-cache-dir
  sudo pip2 install flake8 --upgrade --no-cache-dir
else
  sudo ldconfig
  sudo apt-get install python-dev python-pip python-all-dev
  sudo apt-get install python-setuptools python-wheel
  sudo apt-get install swig
  sudo pip install numpy --upgrade
  sudo pip install scipy --upgrade
  sudo pip install mpi4py --upgrade
  sudo pip install flake8 --upgrade
fi

test -n $CC  && unset CC
test -n $FCC  && unset FCC
test -n $CXX && unset CXX
