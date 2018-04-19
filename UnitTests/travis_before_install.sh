
if [[ "$TESTOS" == "LINUX" ]]; then
  sudo apt-get install gfortran
fi
bash UnitTests/travis_build_mpi.sh
cd openmpi-3.0.1
sudo make install >/dev/null 2>&1
cd ../
if [[ "$TESTOS" == "OSX" ]]; then
  python install python3
  brew install doxygen
  brew install cmake
  brew install pip
else
  sudo ldconfig
  sudo apt-get install python-dev
  sudo apt-get install doxygen
fi
wget https://downloads.sourceforge.net/swig/swig-3.0.12.tar.gz
tar xvf swig-3.0.12.tar.gz >/dev/null
cd swig-3.0.12
./configure >/dev/null 2>&1
make >/dev/null 2>&1
sudo make install >/dev/null 2>&1
cd ../
sudo pip install scipy --upgrade
sudo pip install numpy --upgrade
sudo pip install mpi4py --upgrade
test -n $CC  && unset CC
test -n $FCC  && unset FCC
test -n $CXX && unset CXX
