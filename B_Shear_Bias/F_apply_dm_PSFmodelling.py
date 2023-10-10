# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-01-22 15:11:12
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-22 17:27:34

### apply the dm from PSF modelling to the whole sample results

import os
import re
import glob
import argparse

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

# +++++++++++++++++++++++++++++ general info

# number of bins
N_Zbins = 5

# the m results after varCorr
inpath_m = './results/m_skills_v07D7p1_LF_321_kidsPhotometry_A12_gold_Rewei_varCorr_dz0p1_2D.csv'

# the dm results from PSF modelling (without gold selection)
inpath_dm = '/net/grecht/data2/ssli_files/Projects/Projects/8ShearBias_ImSim/MultiBand_ImSim/correction_varShear_PSFmodelling\
/results/dm_PSFmodelling_41_rewei.csv'
# where to save
outpath = inpath_m.replace('.csv', '_PSFmodellingCorr.csv')

# # the dm results from PSF modelling (with gold selection)
# inpath_dm = './results/dm_PSFmodelling_41_rewei_gold.csv'
# # where to save
# outpath = inpath_m.replace('.csv', '_PSFmodellingCorr_gold.csv')

# +++++++++++++++++++++++++++++ load catalogues
# load m results
cata_m = pd.read_csv(inpath_m)[:(N_Zbins+1)]
cata_dm = pd.read_csv(inpath_dm)[:(N_Zbins+1)]

# +++++++++++++++++++++++++++ correct the m 
for col_m in ['m1', 'm2']:
    cata_m.loc[:, col_m] += cata_dm.loc[:, col_m]
    cata_m.loc[:, f'{col_m}_err'] = (cata_m.loc[:, f'{col_m}_err']**2 
                                    + cata_dm.loc[:, f'{col_m}_err']**2
                                    )**0.5  
# save
cata_m.to_csv(outpath, index=False)
print(cata_m)
print('saved to', outpath)