set -e

if [[ "$TESTOS" == "LINUX" ]]; then
   conda activate ntpoly-conda
   which python
   conda env list
fi

cd Build
export CTEST_OUTPUT_ON_FAILURE=1
eval "$MAKETEST"
cd ../
