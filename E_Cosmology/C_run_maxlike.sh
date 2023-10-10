#!/usr/bin/bash

# @Author: lshuns
# @Date:   2023-05-01 10:39:51
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-21 17:00:54

### run local minimisation with Nelder-Mead method

# activate kcap
source /disks/shear10/ssli/K1000CS/my_kcap/env/bin/activate
source cosmosis-configure

# how many cores for each run
export OMP_NUM_THREADS=100

## loop over all runs
# TAGs=("multinest_nIA" "multinest_wIA" "multinest_wIA_eta" "multinest_wIA_noCut")
#### hydra: Elapsed:6:22:41.31,User=803337.805,System=1120.597,CPU=3503.5%.
# TAGs=("multinest_nIA_noCut")
TAGs=("polychord_nIA")

for TAG in "${TAGs[@]}"
do
    # polychord and multinest have different prefix
    if [[ ${TAG} == *"polychord"* ]]; then
        PREFIX="chain"
    else
        PREFIX="output"
    fi

    # run Marika's code
    python ./maxlike/maxlike_cosmosis.py -i ./LF321_${TAG}/pipeline.ini \
      -m /disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/${TAG}/${PREFIX}_${TAG}.txt\
      -o /disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/${TAG}/maxpost_${TAG}_start.txt\
      --max_post -s \
      --maxiter 3000
done
