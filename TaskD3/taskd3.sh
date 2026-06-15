#!/usr/bin/env bash
set -euo pipefail

tar -xzf myenv.tar.gz
source myenv/bin/activate
mkdir -p saved-calcs_d3
#python -m pip install --quiet -r requirements.txt
which python
# Set the license path and home directory.
HOME=$_CONDOR_SCRATCH_DIR
#python /arc-basis.py $1
export MOSEKLM_LICENSE_FILE="$PWD/mosek/mosek.lic"
export HOME="$_CONDOR_SCRATCH_DIR"

myenv/bin/python arc-basis_d3.py "$1" "$2"

#tar -czf saved_calcs/saved-calcs_$1_$2.tar.gz saved_calcs/saved-calcs_$1_$2
#tar -czf data$1.tar.gz saved_data/