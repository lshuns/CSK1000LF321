# -*- coding: utf-8 -*-
# @Author: ssli
# @Date:   2022-10-27 21:24:13
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-22 22:11:39

### split the catalogue to unique shear values

import os
import re
import glob
import argparse

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

# +++++++++++++++++++++++++++++ general info

# input 
inpath = '/disks/shear10/ssli/K1000CS/skills_v07D7p1_Outputs/skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_newCut_A1_WeiCut_goldclasses.cat.A2.feather'

# shear and rotation info
unique_shear_tags = ['m283m283', 'm283p283', 'p283m283', 'p283p283']
# unique_rots = [0, 90]

# +++++++++++++++++++++++++++++ workhorse

# load the catalogue
cata_final = pd.read_feather(inpath)
print('number in ori', len(cata_final))

# select shear
Ntot = 0
for run_tag in unique_shear_tags:

    if 'unique_rots' in locals():
        for rot in unique_rots:
            mask_shear = cata_final['run_tag'].values == run_tag
            mask_rot = cata_final['gal_rot'].values == rot

            cata_selec = cata_final[mask_shear&mask_rot]
            del mask_shear, mask_rot

            cata_selec.reset_index(drop=True, inplace=True)
            print(f'number with shear {run_tag}, {rot}:', len(cata_selec))
            Ntot += len(cata_selec)

            # save
            out_path = inpath + f'.shear_{run_tag}_rot_{rot:.0f}.feather'
            ## check existence
            if os.path.isfile(out_path):
                os.remove(out_path)

            cata_selec.to_feather(out_path)
            print(f'catalogue saved as {out_path}')
            del cata_selec
    else:
        mask_shear = cata_final['run_tag'].values == run_tag
        cata_selec = cata_final[mask_shear]
        del mask_shear
        cata_selec.reset_index(drop=True, inplace=True)
        print(f'number with shear {run_tag}:', len(cata_selec))
        Ntot += len(cata_selec)
        # save
        out_path = inpath + f'.shear_{run_tag}.feather'
        ## check existence
        if os.path.isfile(out_path):
            os.remove(out_path)
        cata_selec.to_feather(out_path)
        print(f'catalogue saved as {out_path}')
        del cata_selec

print('number in split', Ntot)

# number in ori 28885361
# number with shear m283m283: 7219165
# catalogue saved as /disks/shear10/ssli/K1000CS/skills_v07D7p1_Outputs/skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_newCut_A1_WeiCut_goldclasses.cat.A2.feather.shear_m283m283.feather
# number with shear m283p283: 7222023
# catalogue saved as /disks/shear10/ssli/K1000CS/skills_v07D7p1_Outputs/skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_newCut_A1_WeiCut_goldclasses.cat.A2.feather.shear_m283p283.feather
# number with shear p283m283: 7221808
# catalogue saved as /disks/shear10/ssli/K1000CS/skills_v07D7p1_Outputs/skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_newCut_A1_WeiCut_goldclasses.cat.A2.feather.shear_p283m283.feather
# number with shear p283p283: 7222365
# catalogue saved as /disks/shear10/ssli/K1000CS/skills_v07D7p1_Outputs/skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_newCut_A1_WeiCut_goldclasses.cat.A2.feather.shear_p283p283.feather
# number in split 28885361
# Elapsed:1:14.85,User=213.305,System=127.636,CPU=455.4%.
