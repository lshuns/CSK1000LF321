# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-10-17 09:52:10
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-12 13:46:20

### apply photometry & shear selections
###### Reference: Appendix D of Li et al. (2022)
###### without weight cut

import os
import re
import glob
import argparse

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

# +++++++++++++++++++++++++++++ general info

# in info
inpath = '/disks/shear10/ssli/KiDS/K1000_LF_321_mosaic/KIDS_1006tiles_r_SDSS.V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2.feather'
# out info
outpath = inpath.replace('.feather', '_selec_noWeiCut.feather')

# +++++++++++++++++++++++++++++ workhorse

## load and select catalogue
cata = pd.read_feather(inpath)
print('number original', len(cata))

# >>>>>>>>>> mask
nine_band_mask = (cata['MASK'].values & 28668)==0

# >>>>>>>>>> photometry-related
## 1. 9-band photometry cut
flag_9 = np.zeros(len(cata)) 
for band in ['u', 'g', 'r', 'i', 'Z', 'Y', 'J', 'H', 'Ks']:
    flag_9 += cata[f'FLAG_GAAP_{band}'].values
mask_gaap = (flag_9==0)
del flag_9

## 2. remove asteroids
gmr = np.array(cata['MAG_GAAP_g']) - np.array(cata['MAG_GAAP_r'])
imr = np.array(cata['MAG_GAAP_i']) - np.array(cata['MAG_GAAP_r'])
mask_ast = (gmr <= 1.5) | (imr <= 1.5)
del gmr, imr

# >>>>>>>>>> LF-related
#### a) remove unmeasured
mask_psf = (cata['PSF_Q11'] != 0.0) & (cata['PSF_Q22'] != 0.0)

#### b) remove binaries
mask_binary = (np.hypot(cata['autocal_e1'].values, cata['autocal_e2'].values) <= 0.8) \
            | (cata['autocal_scalelength_pixels'] >= \
                (0.5 * np.exp(0.65788*(24.2 - cata['MAG_GAAP_r']))))

#### c) fitclass cut
mask_class = (cata['fitclass']!=-1) \
            & (cata['fitclass']!=-10) \
            & (cata['fitclass']!=-4) \
            & (cata['fitclass']!=1) \
            & (cata['fitclass']!=2) \
            & (cata['fitclass']!=-7) \
            & (cata['fitclass']!=-3)

#### d) magnitude cut
mask_mag = (cata['MAG_AUTO']>20.0)

#### e) blending cut
mask_blending = (cata['contamination_radius']>4.25)

## apply
cata = cata[nine_band_mask & mask_gaap & mask_ast &
            mask_psf & mask_binary & mask_class & mask_mag & mask_blending]
cata.reset_index(drop=True, inplace=True)
del nine_band_mask
del mask_gaap, mask_ast
del mask_psf, mask_binary, mask_class, mask_mag, mask_blending
print('>>> number after fiducial selection', len(cata))

# >>>>>>>>>> new selection
## 1. avoid negative variance or snr
cata = cata[(cata['2D_measurement_variance'].values>0)&(cata['model_SNratio'].values>0)]
cata.reset_index(drop=True, inplace=True)

## 2. resolution cut
### circularised galaxy size
emod = np.hypot(cata['autocal_e1'].values, cata['autocal_e2'].values)
cata.loc[:, 'r_ab'] = np.array(cata['autocal_scalelength_pixels'].values) * np.sqrt((1.-emod)/(1.+emod))
del emod
### PSF size
cata.loc[:, 'PSFsize'] = np.array(
                (cata['PSF_Q11'].values*cata['PSF_Q22'].values \
                    - cata['PSF_Q12'].values**2.)**0.5)
### resolution parameter
cata.loc[:, 'R'] = np.array(
                cata['PSFsize'].values\
                    / (cata['r_ab'].values**2 + cata['PSFsize'].values))
### apply
cata = cata[cata['R']<0.9]
cata.reset_index(drop=True, inplace=True)

## 3. size cut
cata = cata[((cata['raw_scalelength_pixels'].values)>=0.5)]
cata.reset_index(drop=True, inplace=True)

print('>>> number after new selection', len(cata))

# >>>>>>>>>> save
cata.to_feather(outpath)
print('final cata saved to', outpath)

# number original 112,907,874
# >>> number after fiducial selection 48,428,603
# >>> number after new selection 38,026,162
# final cata saved to /disks/shear10/ssli/KiDS/K1000_LF_321_mosaic/KIDS_1006tiles_r_SDSS.V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_selec_noWeiCut.feather
# Elapsed:5:44.54,User=371.924,System=1281.973,CPU=480.0%.
