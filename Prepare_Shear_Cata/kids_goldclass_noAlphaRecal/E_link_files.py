# -*- coding: utf-8 -*-
# @Author: ssli
# @Date:   2022-11-02 15:47:21
# @Last Modified by:   ssli
# @Last Modified time: 2022-12-24 12:27:52

### link the file for CosmoWrapper

import os 

# 1. the specz sample
inpath = '/disks/shear10/ssli/K1000CS/K1000_info/KiDS_specz_PAUS_COSMOS2015.csv'
outdir = '/disks/shear10/ssli/K1000CS/LF321_noAlphaRecal_Inputs'
os.symlink(inpath, os.path.join(outdir, os.path.basename(inpath)))
print('symlink created for', inpath)

# 2. the shear catalogue
inpath = '/disks/shear16/ssli/KiDS/K1000_LF_321_mosaic/KIDS_1006tiles_r_SDSS.V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_selec_noWeiCut_WeiCut.cat'
outpath = '/disks/shear10/ssli/K1000CS/LF321_noAlphaRecal_Inputs/K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_noA.cat'
os.symlink(inpath, outpath)
print('symlink created as', outpath)
