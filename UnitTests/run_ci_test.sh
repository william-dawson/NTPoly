set -e

if [[ "$TESTOS" == "LINUX" ]]; then
   conda activate ntpoly-conda
fi

if [[ "$THREADOFF" == "1" ]]; then
   export OMP_NUM_THREADS=1
   echo "Setting threads to 1"
else
   export OMP_NUM_THREADS=2
   echo "Setting threads to 2"
fi

cd Build
export CTEST_OUTPUT_ON_FAILURE=1
eval "$MAKETEST"
cd ../
