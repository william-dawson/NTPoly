
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
elif [[ "$TESTOS" == "OSX" ]]; then
  brew reinstall gcc
  brew link --overwrite gcc
  brew install open-mpi
  brew install doxygen
  brew install cmake
  brew install clang-format
  brew install emacs
  brew install ninja
fi

conda activate
conda env create -f environment.yml
