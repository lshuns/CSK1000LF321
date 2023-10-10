#!/usr/bin/bash

# @Author: lshuns
# @Date:   2023-04-27 21:21:45
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-12 09:44:34

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


# ######## markermerr
# ...
# sample_S8 took: 0.000 seconds
# sigma8toAs took: 1.919 seconds
# camb took: 46.179 seconds
# extrapolate took: 0.097 seconds
# correlated_dz_priors took: 0.000 seconds
# fits_nz took: 0.000 seconds
# photoz_bias took: 0.004 seconds
# linear_alignment took: 0.003 seconds
# pk_to_cl took: 2.341 seconds
# add_intrinsic took: 0.004 seconds
# cosebis took: 1.953 seconds
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
# Total pipeline time: 52.5 seconds
# Pipeline ran okay.
#     Likelihood cosebis = -29.464036400982486
# Likelihood total = -29.464036400982486
# Saving 13642 samples
# Saving 13642 samples
# Elapsed:13:48:24.28,User=2979199.178,System=2606.210,CPU=5999.0%.
