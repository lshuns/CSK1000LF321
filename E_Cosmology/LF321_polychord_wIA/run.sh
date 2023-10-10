#!/usr/bin/bash

# @Author: lshuns
# @Date:   2023-04-27 21:21:45
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-10 22:13:39

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
