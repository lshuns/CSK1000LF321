# -*- coding: utf-8 -*-
# @Author: ssli
# @Date:   2022-10-24 09:24:59
# @Last Modified by:   ssli
# @Last Modified time: 2022-12-24 15:01:02

### link the file for SOM gold selection

import os 

# 1. the specz sample
inpath = '/disks/shear10/ssli/K1000CS/K1000_info/KiDS_specz_PAUS_COSMOS2015.csv'
outdir = '/disks/shear10/ssli/K1000CS/skills_v07D7_Inputs'
os.symlink(inpath, os.path.join(outdir, os.path.basename(inpath)))
print('symlink created for', inpath)