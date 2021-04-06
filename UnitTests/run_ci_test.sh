set -e

if [[ "$CONDA" == "1" ]]; then
   conda activate ntpoly-conda
fi

cd Build
export CTEST_OUTPUT_ON_FAILURE=1
eval "$MAKETEST"
cd ../
