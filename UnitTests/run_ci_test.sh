set -e

if [[ "$TESTOS" == "LINUX" ]]; then
   conda activate ntpoly-conda
fi

echo "Setting threads to $OMP_NUM_THREADS"

cd Build
export CTEST_OUTPUT_ON_FAILURE=1
eval "$MAKETEST"
cd ../
