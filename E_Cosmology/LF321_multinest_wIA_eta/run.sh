#!/usr/bin/bash

# @Author: lshuns
# @Date:   2023-04-27 21:21:45
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-13 15:40:56

# for mpi run
module load OpenMPI

# activate kcap
source /disks/shear10/ssli/K1000CS/my_kcap/env/bin/activate
source cosmosis-configure

# how many runs
NUM_PROCESSES=100

# how many cores for each run
export OMP_NUM_THREADS=1

# start run
mpirun -n ${NUM_PROCESSES} cosmosis --mpi pipeline.ini

########### hydra
# sample_S8 took: 0.000 seconds
# sigma8toAs took: 1.210 seconds
# camb took: 28.993 seconds
# extrapolate took: 0.057 seconds
# correlated_dz_priors took: 0.000 seconds
# fits_nz took: 0.000 seconds
# photoz_bias took: 0.003 seconds
# linear_alignment took: 0.002 seconds
# pk_to_cl took: 1.352 seconds
# add_intrinsic took: 0.003 seconds
# cosebis took: 1.412 seconds
# Gathering theory outputs to make a vector
#   Skipped xip
#   Skipped xim
#   Skipped PneE
#   Skipped PeeE
#   cosebis stops at bin 5
# Masking 0 values in En because they had ell or theta outside (0.5,5.5)
#   Did scale cuts to theory
# scale_cuts took: 0.001 seconds
# cosebis_like took: 0.000 seconds
# Total pipeline time: 33.0 seconds
# Pipeline ran okay.
#     Likelihood cosebis = -30.079252683436295
# Likelihood total = -30.079252683436295
# Saving 14928 samples
# Saving 14928 samples
# Elapsed:5:32:37.80,User=1993017.816,System=1097.835,CPU=9991.6%.
