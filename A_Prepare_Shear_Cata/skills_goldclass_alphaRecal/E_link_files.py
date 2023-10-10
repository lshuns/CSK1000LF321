# -*- coding: utf-8 -*-
# @Author: ssli
# @Date:   2022-10-24 09:24:59
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-22 12:36:32

### link the file for SOM gold selection

import os 

# 1. the adapted specz sample
inpath = '/disks/shear10/ssli/K1000CS/K1000_info/KiDS_specz_PAUS_COSMOS2015_adapt.fits'
outdir = '/disks/shear10/ssli/K1000CS/skills_v07D7p1_Outputs'
os.symlink(inpath, os.path.join(outdir, os.path.basename(inpath)))
print('symlink created for', inpath)

# 3. the SOM
inpath = '/disks/shear10/ssli/K1000CS/K1000_info/Spec_Train_Adapt_DIR_SOMdata.Rdata'
outdir = '/disks/shear10/ssli/K1000CS/skills_v07D7p1_Outputs'
os.symlink(inpath, os.path.join(outdir, os.path.basename(inpath)))
print('symlink created for', inpath)