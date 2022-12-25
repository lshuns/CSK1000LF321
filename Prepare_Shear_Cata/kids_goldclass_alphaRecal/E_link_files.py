# -*- coding: utf-8 -*-
# @Author: ssli
# @Date:   2022-09-24 11:33:32
# @Last Modified by:   ssli
# @Last Modified time: 2022-12-24 11:26:06

### link the file for SOM gold selection

import os 

# 1. the specz sample
inpath = '/disks/shear10/ssli/K1000CS/K1000_info/KiDS_specz_PAUS_COSMOS2015.csv'
outdir = '/disks/shear10/ssli/K1000CS/LF321_Inputs'
os.symlink(inpath, os.path.join(outdir, os.path.basename(inpath)))
print('symlink created for', inpath)

# 2. the shear catalogue
inpath = '/disks/shear16/ssli/KiDS/K1000_LF_321_mosaic/KIDS_1006tiles_r_SDSS.V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_selec_A1_WeiCut.cat'
outpath = '/disks/shear10/ssli/K1000CS/LF321_Inputs/K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A1.cat'
os.symlink(inpath, outpath)
print('symlink created as', outpath)
