set -e

if [[ "$CONDA" == "1" ]]; then
   ls ~
   echo $CONDA
   source $CONDA/etc/profile.d/conda.sh
   conda activate ntpoly-conda-env
fi

cd Build
export CTEST_OUTPUT_ON_FAILURE=1
eval "$MAKETEST"
cd ../
