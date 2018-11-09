
if [[ "$TESTOS" == "LINUX" ]]; then
  sudo apt-get install gfortran
  sudo apt-get install libblas-dev liblapack-dev
  sudo apt-get install openmpi-bin openmpi-doc libopenmpi-dev
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
else
  sudo ldconfig
  sudo apt-get install python-dev
  if [[ -z ${SKIPDOXYGEN+x} ]]; then
    sudo apt-get install doxygen
  fi
  if [[ -z ${SKIPSWIG+x} ]] ; then
    if [ -f "swig-3.0.12/README" ]; then
      echo "Using cached swig";
    else
      wget https://downloads.sourceforge.net/swig/swig-3.0.12.tar.gz;
      tar xvf swig-3.0.12.tar.gz >/dev/null;
    fi
    cd swig-3.0.12;
    ./configure >/dev/null 2>&1;
    make >/dev/null 2>&1;
    sudo make install >/dev/null 2>&1
    cd ../
  fi
  sudo pip install scipy --upgrade
  sudo pip install numpy --upgrade
  sudo pip install mpi4py --upgrade
fi

test -n $CC  && unset CC
test -n $FCC  && unset FCC
test -n $CXX && unset CXX
