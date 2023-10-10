#!/usr/bin/bash

# @Author: lshuns
# @Date:   2023-05-01 10:39:51
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-21 20:33:05

### get the summary statistics from CosmoSIS

# activate kcap
source /disks/shear10/ssli/K1000CS/my_kcap/env/bin/activate
source cosmosis-configure

# loop over all runs
# TAGs=("multinest_K1000_like" "multinest_nIA" "multinest_nIA_noCut" "multinest_wIA" "multinest_wIA_eta" "multinest_wIA_noCut")
TAGs=("polychord_nIA")

for TAG in "${TAGs[@]}"
do

    # polychord and multinest have different prefix
    if [[ ${TAG} == *"polychord"* ]]; then
        PREFIX="chain"
    else
        PREFIX="output"
    fi

    # run 
    cosmosis-postprocess --no-plots -o /disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/${TAG}/cosmosis_summary/ \
        /disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/${TAG}/${PREFIX}_${TAG}.txt
done
