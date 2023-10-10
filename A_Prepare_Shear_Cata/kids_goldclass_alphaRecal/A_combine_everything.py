# -*- coding: utf-8 -*-
# @Author: ssli
# @Date:   2022-10-17 09:46:17
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-12 13:45:09

### combine K1000-321 superuser catalogues
###### everything without selections
###### select columns for memory

import os
import glob

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

######################## I/O

# csv file includes wanted columns
wanted_cols = '../utils/wanted_cols_LF321.csv'

# superuser catalogues
inpath_list = glob.glob(os.path.join('/disks/shear10/ssli/KiDS/K1000_LF_321', 'KIDS_*_r_SDSS.V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2.cat'))
Ntiles = len(inpath_list)
print('number of tiles found', Ntiles)

# out files
outpath = f'/disks/shear10/ssli/KiDS/K1000_LF_321_mosaic/KIDS_{Ntiles}tiles_r_SDSS.V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2.feather'

######################## work horse

# 1. wanted info
cols_final = pd.read_csv(wanted_cols)['col'].to_list()

# 2. get catalogues
cata_final = []
for inpath in inpath_list:

    # load
    with fits.open(inpath) as hdul:
        cata = Table(hdul[1].data).to_pandas()
    cata = cata[cols_final]

    # collect
    cata_final.append(cata)
    del cata

cata_final = pd.concat(cata_final, ignore_index=True)
print(f'Total number of source {len(cata_final)}')
cata_final.to_feather(outpath)
print('combined cata saved to', outpath)

# number of tiles found 1006
# Total number of source 112,907,874
# combined cata saved to /disks/shear10/ssli/KiDS/K1000_LF_321_mosaic/KIDS_1006tiles_r_SDSS.V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2.feather
# Elapsed:36:28.70,User=1840.703,System=1599.687,CPU=157.1%.
