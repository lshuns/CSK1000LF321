#!/usr/bin/bash

# @Author: lshuns
# @Date:   2023-05-01 10:39:51
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-21 17:05:56

### run local minimisation with Nelder-Mead method
####### to find the S8 shift due to the shift in data vector

# activate kcap
source /disks/shear10/ssli/K1000CS/my_kcap/env/bin/activate
source cosmosis-configure

# how many cores for each run
# export OMP_NUM_THREADS=100

# RUNNAME="multinest_nIA"
RUNNAME="polychord_nIA"
variations=("no_m_bias" "with_raw_m_bias" "test_m_nU" "test_m_nD" "test_m_qU" "test_m_qD" "test_m_sizeU" "test_m_sizeD")

#### multinest
for variation in "${variations[@]}"; do
    # run Marika's code
    python ../maxlike/maxlike_cosmosis.py -i ./config/${RUNNAME}_pipeline_${variation}.ini \
      -m /disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/${RUNNAME}/maxpost_${RUNNAME}_start.txt\
      -o /disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/${RUNNAME}/maxpost_${RUNNAME}_${variation}.txt\
      --max_post -s \
      --maxiter 3000
done