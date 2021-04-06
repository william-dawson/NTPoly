set -e

if [[ "$CONDA" == "1" ]]; then
   source ~/.bashrc
   conda activate ntpoly-conda-env
fi

cd Build
export CTEST_OUTPUT_ON_FAILURE=1
eval "$MAKETEST"
cd ../
