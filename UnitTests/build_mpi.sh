# FOR TRAVIS-CI ONLY
if [ -f "openmpi-2.1.1/README" ]; then
  echo "Using cached openmpi";
else
  wget https://www.open-mpi.org/software/ompi/v2.1/downloads/openmpi-2.1.1.tar.gz;
  tar xvf openmpi-2.1.1.tar.gz >/dev/null;
  cd openmpi-2.1.1;
  ./configure >/dev/null 2>&1;
  make >/dev/null 2>&1;
fi
