# -*- coding: utf-8 -*-
# @Author: ssli
# @Date:   2022-11-02 15:26:26
# @Last Modified by:   ssli
# @Last Modified time: 2022-12-24 15:46:54

### apply weight cut and transfer to LDAC cat
###### for catalogues without alphaRecal

import os

import pandas as pd

from astropy.io import fits
from astropy.table import Table

# +++++++++++++ I/O

# the original catalogue
inpath = '/disks/shear16/ssli/KiDS/K1000_LF_321_mosaic/KIDS_1006tiles_r_SDSS.V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_selec_noWeiCut.feather'

# weight column
wei_col = 'weight'

# out name
outpath = os.path.join(os.path.dirname(inpath), os.path.basename(inpath).replace('_noWeiCut_', '_').replace('.feather', '_WeiCut.cat'))

# the LDAC template
ldac_template = "../../utils/LDAC.cat"

# csv file includes wanted columns
wanted_cols = '../utils/wanted_cols_LF321_noA.csv'

# +++++++++++++ workhorse

# load the data file
cata = pd.read_feather(inpath)
print('number original', len(cata))

# weight cut
cata = cata[cata[wei_col]>0]
cata.reset_index(drop=True, inplace=True)
print('number useful', len(cata))

# to hdu
## column names and format
tmp = pd.read_csv(wanted_cols)
col_names = tmp['col'].to_list()
col_format = tmp['format'].to_list()
del tmp
## build the col in hdul format
cols_list = []
for i_col, col in enumerate(col_names):
    cols_list.append(fits.Column(name=col, array=cata[col].values, format=col_format[i_col]))
del cata
hdu = fits.BinTableHDU.from_columns(cols_list)
del cols_list
hdu.name = "OBJECTS"
## add dummy TCOMM values
for i, col in enumerate(col_names):
    hdu.header.set('TCOMM%d' % (i + 1), col,
                           after='TTYPE%d' % (i + 1))
del col_names, col_format

# combine with template LDAC HDU
with fits.open(ldac_template) as ldac:
    ldac.append(hdu)
    del hdu
    print("LDAC info:")
    ldac.info()
    ldac.writeto(outpath, overwrite=True)
print(f"saved as {outpath}")

# number original 38,026,162
# number useful 32,180,207
# LDAC info:
# Filename: ../../utils/LDAC.cat
# No.    Name      Ver    Type      Cards   Dimensions   Format
#   0  PRIMARY       1 PrimaryHDU       4   ()      
#   1  FIELDS        1 BinTableHDU     24   65534R x 5C   [1J, 1J, 1J, 16A, 1J]   
#   2  OBJECTS       1 BinTableHDU    342   32180207R x 111C   ['16A', '1J', '1J', '1E', '1E', '1D', '1D', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1J', '1E', '1I', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E']   
# saved as /disks/shear16/ssli/KiDS/K1000_LF_321_mosaic/KIDS_1006tiles_r_SDSS.V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_selec_noWeiCut_WeiCut.cat
# Elapsed:9:31.41,User=435.141,System=531.457,CPU=169.1%.
