# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-04-18 16:29:58
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-24 10:55:47

### calculate the neff & ellipticity dispersion

import os

import numpy as np
import pandas as pd

from astropy.io import fits

# >>>>>>>>>>> I/O

# the input catalogue
inpath = '/disks/shear10/ssli/K1000CS/LF321_Outputs/K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A1_goldclasses.cat.A2.feather'
cata_type = 'feather' 
col_e1 = 'AlphaRecalD2_e1'
col_e2 = 'AlphaRecalD2_e2'
col_wei = 'AlphaRecalC_weight'
col_ZB = 'Z_B'
col_SOM = 'Flag_SOM_Fid_NONE'

# where to save
outdir = '/disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/COV/input'

# K1000 effective area from healpix
SurveyArea = 3.1212e+06 # 867.0 deg2

# K1000 tomographic binning
ZB_ranges = (0.1, 0.3, 0.5, 0.7, 0.9, 1.2)

# >>>>>>>>>>>> workhorse

# load cata
if cata_type == 'ldac':
    with fits.open(inpath) as hdul:
        cata = hdul['OBJECTS'].data
    print('number ori', len(cata))
    # we only want gold classes
    cata = cata[(cata[col_SOM]==1)&(cata[col_wei]>0)]
    print('number selected', len(cata))
    # used values
    e1 = cata[col_e1]
    e2 = cata[col_e2]
    weight = cata[col_wei]
    Z_B = cata[col_ZB]
    del cata
elif cata_type == 'feather':
    cata = pd.read_feather(inpath)
    print('number ori', len(cata))
    # we only want gold classes and non-zero weight
    cata = cata[(cata[col_SOM]==1)&(cata[col_wei]>0)]
    print('number selected', len(cata))
    # used values
    e1 = cata[col_e1]
    e2 = cata[col_e2]
    weight = cata[col_wei]
    Z_B = cata[col_ZB]
    del cata

# loop over tomo bins
n_eff_list = np.zeros(len(ZB_ranges) - 1)
sigma_e_list = np.zeros(len(ZB_ranges) - 1)
for i_zbin in range(len(ZB_ranges) - 1):

    # 0.001 is K1000 convention, have no idea why
    low_z_cut = ZB_ranges[i_zbin] + 0.001
    high_z_cut = ZB_ranges[i_zbin + 1] + 0.001

    # select 
    mask_tmp = (Z_B > low_z_cut) & (Z_B <= high_z_cut)
    weight_masked = weight[mask_tmp]
    weight_sq_masked = np.square(weight_masked)
    e1_masked = e1[mask_tmp]
    e2_masked = e2[mask_tmp]
    del mask_tmp

    # effective number density
    n_eff_list[i_zbin] = np.sum(weight_masked)**2 / np.sum(weight_sq_masked) / SurveyArea

    # e dispersion
    c1 = np.average(e1_masked, weights=weight_masked)
    c2 = np.average(e2_masked, weights=weight_masked)
    sigma_e_list[i_zbin] = (  # geometric mean
        np.sqrt(  # standard deviation of e1
          np.sum(np.square(weight_masked * (e1_masked - c1))) /
          np.sum(weight_sq_masked)
        ) / 2.0 +
        np.sqrt(  # standard deviation of e2
          np.sum(np.square(weight_masked * (e2_masked - c2))) /
          np.sum(weight_sq_masked)
        ) / 2.0)

    del weight_masked, weight_sq_masked, e1_masked, e2_masked

# save the results
outpath = os.path.join(outdir, 'sigmae.txt')
np.savetxt(outpath, sigma_e_list, fmt='%.12f', newline=" ")
print('>>> sigma_e', sigma_e_list)
print('saved to', outpath)
outpath = os.path.join(outdir, 'neff.txt')
np.savetxt(outpath, n_eff_list, fmt='%.12f', newline=" ")
print('>>> n_eff', n_eff_list)
print('saved to', outpath)