# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-04-18 19:50:00
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-03 11:30:11

### link to the XI measurements with the naming used by covariance calculation

import os 

# where to find the XI measurements
indir = '/disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/'
# theta_tag = 'theta_0.5_300.0'
theta_tag = 'theta_2.0_300.0'
in_prefix = f'XI_K1000_ALL_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_{theta_tag}_zbins'

# where to save
outdir = '/disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/COV/input/'
out_prefix = f'XI_K1000_ALL_BLIND_NONE_{theta_tag}'

# number of tomo bins
N_Zbins = 5

# loop over all and link
for i_zbin in range(N_Zbins):
    for j_zbin in range(i_zbin, N_Zbins):
        inpath = os.path.join(indir, f'{in_prefix}_{i_zbin+1}_{j_zbin+1}.asc')
        if not os.path.exists(inpath):
            raise Exception(f'{inpath} does not exist!')

        outpath = os.path.join(outdir, f'{out_prefix}_nBins_5_Bin{i_zbin+1}_Bin{j_zbin+1}.asc')

        os.symlink(inpath, outpath)
        print('symlink created to', outpath)