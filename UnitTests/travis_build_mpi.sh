## @file travis_build_mpi.sh
## @brief FOR TRAVIS-CI ONLY
## This downloads and builds openmpi from the website. This is necessary
## since Travis-CI version of linux is too old to have version 2+ of openmpi.
## This script checks for a cached version before downloading.
if [ -f "openmpi-3.0.1/README" ]; then
  echo "Using cached openmpi";
else
  wget https://www.open-mpi.org/software/ompi/v3.0/downloads/openmpi-3.0.1.tar.gz;
  tar xvf openmpi-3.0.1.tar.gz >/dev/null;
  cd openmpi-2.1.1;
  ./configure >/dev/null 2>&1;
  make >/dev/null 2>&1;
fi
