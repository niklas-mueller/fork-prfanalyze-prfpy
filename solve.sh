#! /bin/bash

# all we have to do is exec python...
export PRF_SOLVER="prfpy"
. /opt/conda/etc/profile.d/conda.sh
conda activate prfpy_analysis
echo "-------------- Activated conda environment ---------------"
exec python /scripts/run_prfpy.py "$1" "$2" "$3" "$4" "$5"