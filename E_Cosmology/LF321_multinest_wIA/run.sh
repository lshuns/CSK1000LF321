#!/usr/bin/bash

# @Author: lshuns
# @Date:   2023-04-27 21:21:45
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-13 15:41:44

# for mpi run
module load OpenMPI

# activate kcap
source /disks/shear10/ssli/K1000CS/my_kcap/env/bin/activate
source cosmosis-configure

# how many runs
NUM_PROCESSES=60

# how many cores for each run
export OMP_NUM_THREADS=1

# start run
mpirun -n ${NUM_PROCESSES} cosmosis --mpi pipeline.ini

############## markermeer
# sample_S8 took: 0.000 seconds
# sigma8toAs took: 2.346 seconds
# camb took: 56.047 seconds
# extrapolate took: 0.130 seconds
# correlated_dz_priors took: 0.000 seconds
# fits_nz took: 0.000 seconds
# photoz_bias took: 0.006 seconds
# linear_alignment took: 0.005 seconds
# pk_to_cl took: 2.364 seconds
# add_intrinsic took: 0.006 seconds
# cosebis took: 2.263 seconds
# Gathering theory outputs to make a vector
#   Skipped xip
#   Skipped xim
#   Skipped PneE
#   Skipped PeeE
#   cosebis stops at bin 5
# Masking 0 values in En because they had ell or theta outside (0.5,5.5)
#   Did scale cuts to theory
# scale_cuts took: 0.002 seconds
# cosebis_like took: 0.000 seconds
# Total pipeline time: 63.2 seconds
# Pipeline ran okay.
#     Likelihood cosebis = -30.204230409232455
# Likelihood total = -30.204230409232455
# Saving 15439 samples
# Saving 15439 samples
# Elapsed:16:50:18.90,User=3633107.172,System=3372.215,CPU=5998.9%.
