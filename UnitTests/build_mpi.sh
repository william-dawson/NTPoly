# FOR TRAVIS-CI ONLY
if [ -d "openmpi-2.1.1" ]; then
  echo "Using cached openmpi";
else
  wget https://www.open-mpi.org/software/ompi/v2.1/downloads/openmpi-2.1.1.tar.gz;
  tar xvf openmpi-2.1.1.tar.gz >/dev/null;
fi
