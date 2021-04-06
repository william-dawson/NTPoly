
if [[ "$TESTOS" == "LINUX" ]]; then
  sudo apt-get update
  sudo apt-get install libblas-dev liblapack-dev
  sudo apt-get install gawk
  sudo apt-get install doxygen
  sudo apt-get install clang-format
  sudo apt-get install emacs
  sudo apt-get install ninja-build
  if [[ "$MPICH" == "1" ]]; then
    sudo apt-get install mpich libmpich-dev
  else
    sudo apt-get install openmpi-bin libopenmpi-dev
  fi
fi

if [[ "$TESTOS" == "OSX" ]]; then
  brew reinstall gcc
  brew link --overwrite gcc
  brew install open-mpi
  brew install doxygen
  brew install cmake
  brew install swig
  brew install clang-format
  brew install emacs
  brew install ninja
  sudo pip3 install numpy --upgrade --no-cache-dir
  sudo pip3 install scipy --upgrade --no-cache-dir
  sudo pip3 install mpi4py --upgrade --no-cache-dir
  sudo pip3 install flake8 --upgrade --no-cache-dir
  sudo pip3 install pyyaml --upgrade --no-cache-dir
  sudo pip3 install ford --upgrade --no-cache-dir
else
  sudo apt-get install python-dev python-pip python-all-dev
  sudo apt-get install python-setuptools python-wheel
  sudo apt-get install swig
  sudo pip install numpy --upgrade
  sudo pip install scipy --upgrade
  sudo pip install mpi4py --upgrade
  sudo pip install flake8 --upgrade
  sudo pip install pyyaml --upgrade
  sudo pip install ford --upgrade
fi
