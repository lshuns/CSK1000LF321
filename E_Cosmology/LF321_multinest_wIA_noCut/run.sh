#!/usr/bin/bash

# @Author: lshuns
# @Date:   2023-04-27 21:21:45
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-14 17:03:33

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

######## hydra
# sample_S8 took: 0.000 seconds
# sigma8toAs took: 1.143 seconds
# camb took: 28.108 seconds
# extrapolate took: 0.060 seconds
# correlated_dz_priors took: 0.000 seconds
# fits_nz took: 0.000 seconds
# photoz_bias took: 0.003 seconds
# linear_alignment took: 0.002 seconds
# pk_to_cl took: 1.170 seconds
# add_intrinsic took: 0.003 seconds
# cosebis took: 1.287 seconds
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
# Total pipeline time: 31.8 seconds
# Pipeline ran okay.
#     Likelihood cosebis = -31.047169036435644
# Likelihood total = -31.047169036435644
# Saving 15290 samples
# Saving 15290 samples
# Elapsed:5:33:32.82,User=2198003.911,System=1357.961,CPU=10989.7%.
