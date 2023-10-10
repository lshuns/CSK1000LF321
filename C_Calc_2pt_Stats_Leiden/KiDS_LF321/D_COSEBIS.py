# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-04-15 14:32:32
# @Last Modified by:   lshuns
# @Last Modified time: 2023-06-26 10:20:18

### calculate the COSEBIS
###### This is a hack from https://github.com/KiDS-WL/Cat_to_Obs_K1000_P1/blob/master/Calc_2pt_Stats/doall_calc2pt.sh
######### mode: \"COSEBIS\": calculate En/Bn for tomo bin pair i j "

import os
import subprocess

# >>>>>>>>> I/O

# the XI measurements
indir = "/disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI"
cat_tag = "K1000_ALL_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid"
# suffix_XI = 'nbins_4000_theta_0.5_300.0'
suffix_XI = 'nbins_4000_theta_2.0_300.0'

# the pre-computed COSEBIS tables
SRCLOC = '../cosebis'
## COSEBIS setups
# BP_COSEBIS_THETAINFO = [0.5, 300]
BP_COSEBIS_THETAINFO = [2.0, 300]
nmodes = 20
# nmodes = 5

# where to save
out_dir = f'/disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/COSEBIS/nmodes{nmodes}'

# have we run linear or log binning for the 2pt correlation function?
binning='log'

# the ZB bins 
N_Zbins = 5

# >>>>>>>>>> workhorse

# check that the pre-computed COSEBIS tables exist
normfile = os.path.join(SRCLOC, 
    'TLogsRootsAndNorms', 
    f'Normalization_{BP_COSEBIS_THETAINFO[0]:.2f}-{BP_COSEBIS_THETAINFO[1]:.2f}.table')
if not os.path.exists(normfile):
    raise Exception(f'{normfile} not found!')
rootfile = os.path.join(SRCLOC,
    'TLogsRootsAndNorms',
    f'Root_{BP_COSEBIS_THETAINFO[0]:.2f}-{BP_COSEBIS_THETAINFO[1]:.2f}.table')
if not os.path.exists(rootfile):
    raise Exception(f'{rootfile} not found!')

# check the COSEBIS code is there
COSEBIS_src = os.path.join(SRCLOC, 'run_measure_cosebis_cats2stats.py')
if not os.path.exists(COSEBIS_src):
    raise Exception(f'{COSEBIS_src} not found!')

# loop over all ZB bins
for i_bin1 in range(N_Zbins):
    for i_bin2 in range(i_bin1, N_Zbins):

        tomoPairTag = f'zbins_{i_bin1+1}_{i_bin2+1}'

        print('>>> working on', tomoPairTag)

        # the XI measurements
        treePath = os.path.join(indir, f'XI_{cat_tag}_{suffix_XI}_{tomoPairTag}.asc')
        # check does the correct input xi file exist?
        if not os.path.exists(treePath):
            raise Exception(f'{treePath} not found!')

        # COSEBI output tag
        filetail = f"COSEBIS_{cat_tag}_theta_{BP_COSEBIS_THETAINFO[0]}_{BP_COSEBIS_THETAINFO[1]}_{tomoPairTag}"

        # Now Integrate output from treecorr with COSEBIS filter functions
        # -i = input file
        # -t = treecorr output theta_col - the first column is zero so -t 1 uses the meanR from Treecorr
        # -p = treecorr output xip_col
        # -m = treecorr output xim_col
        # --cfoldername = output directory
        # -o = filename (outputs En_filename.ascii and Bn_filename.ascii)
        # -b = binning "log" or "lin"
        # -n = number of COSEBIS modes
        # -s = COSEBIS minimum theta
        # -l = COSEBIS maximum theta
        # location of the required pre-compution tables
        # --tfoldername = Tplus_minus    # computes/saves the first time it is run for a fixed theta min/max
        # --norm = TLogsRootsAndNorms/Normalization_${tmin}-${tmax}.table
        # --root = TLogsRootsAndNorms/Root_${tmin}-${tmax}.table

        # build command
        cmd = ['python', COSEBIS_src,
                "-i", treePath, "-t", "1", "-p", "3", "-m", "4",
                "--cfoldername", out_dir, "-o", filetail,
                "-b", binning, "-n", str(nmodes), 
                "-s", str(BP_COSEBIS_THETAINFO[0]), "-l", str(BP_COSEBIS_THETAINFO[1]),
                "--tfoldername", os.path.join(out_dir, 'Tplus_minus'),
                "--norm", normfile, "--root", rootfile
                ]

        # run
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        # I am expecting this to have produced two files called
        Enfile = os.path.join(out_dir, f"En_{filetail}.asc")
        Bnfile = os.path.join(out_dir, f"Bn_{filetail}.asc")
        # Did it work?
        if not os.path.exists(Enfile):
            raise Exception(f"{Enfile} was not created!")
        if not os.path.exists(Bnfile):
            raise Exception(f"{Bnfile} was not created!")
        print(f'{Enfile} produced')
        print(f'{Bnfile} produced') 