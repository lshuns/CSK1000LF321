#!/usr/bin/bash

# @Author: lshuns
# @Date:   2023-04-27 21:21:45
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-14 10:47:49

# for mpi run
module load OpenMPI

# activate kcap
source /disks/shear10/ssli/K1000CS/my_kcap/env/bin/activate
source cosmosis-configure

# how many runs
NUM_PROCESSES=110

# how many cores for each run
export OMP_NUM_THREADS=1

# start run
mpirun -n ${NUM_PROCESSES} cosmosis --mpi pipeline.ini

# ####### hydra
# sample_S8 took: 0.000 seconds
# sigma8toAs took: 1.354 seconds
# camb took: 30.454 seconds
# extrapolate took: 0.055 seconds
# correlated_dz_priors took: 0.000 seconds
# fits_nz took: 0.000 seconds
# photoz_bias took: 0.003 seconds
# linear_alignment took: 0.002 seconds
# pk_to_cl took: 1.227 seconds
# add_intrinsic took: 0.003 seconds
# cosebis took: 1.434 seconds
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
# Total pipeline time: 34.5 seconds
# Pipeline ran okay.
#     Likelihood cosebis = -31.207958472575815
# Likelihood total = -31.207958472575815
# Saving 13613 samples
# Saving 13613 samples
# Elapsed:4:31:57.75,User=1792095.761,System=1130.561,CPU=10989.4%.
