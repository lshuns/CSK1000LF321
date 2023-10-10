# -*- coding: utf-8 -*-
# @Author: ssli
# @Date:   2022-09-24 11:33:32
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-22 12:31:57

### link the file for SOM gold selection

import os 

# 1. the shear catalogue
inpath = '/disks/shear10/ssli/KiDS/K1000_LF_321_mosaic/KIDS_1006tiles_r_SDSS.V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_selec_A1_WeiCut.cat'
outpath = '/disks/shear10/ssli/K1000CS/LF321_Inputs/K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A1.cat'
os.symlink(inpath, outpath)
print('symlink created as', outpath)

# 2. the adapted specz sample
inpath = '/disks/shear10/ssli/K1000CS/K1000_info/KiDS_specz_PAUS_COSMOS2015_adapt.fits'
outdir = '/disks/shear10/ssli/K1000CS/LF321_Outputs'
os.symlink(inpath, os.path.join(outdir, os.path.basename(inpath)))
print('symlink created for', inpath)

# 3. the SOM
inpath = '/disks/shear10/ssli/K1000CS/K1000_info/Spec_Train_Adapt_DIR_SOMdata.Rdata'
outdir = '/disks/shear10/ssli/K1000CS/LF321_Outputs'
os.symlink(inpath, os.path.join(outdir, os.path.basename(inpath)))
print('symlink created for', inpath)