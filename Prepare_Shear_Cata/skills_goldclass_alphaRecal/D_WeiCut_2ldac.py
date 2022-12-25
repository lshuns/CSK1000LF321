# -*- coding: utf-8 -*-
# @Author: ssli
# @Date:   2022-10-24 09:24:45
# @Last Modified by:   ssli
# @Last Modified time: 2022-12-24 15:43:55

### apply weight cut and transfer to LDAC cat

import os

import pandas as pd

from astropy.io import fits
from astropy.table import Table

# +++++++++++++ I/O

# the LDAC template
ldac_template = "../../utils/LDAC.cat"

# csv file includes wanted columns
wanted_cols = '../utils/wanted_cols_skills_v07D7_A1.csv'

# the original catalogue
inpath = '/disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_A1.feather'

# weight column
wei_col = 'AlphaRecalC_weight'

# out name
outpath = os.path.join('/disks/shear10/ssli/K1000CS/skills_v07D7_Inputs',
                            os.path.basename(inpath).replace('_noWeiCut_', '_').replace('.feather', '_WeiCut.cat'))

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

# number original 47,870,631
# number useful 39,103,516
# LDAC info:
# Filename: ../../utils/LDAC.cat
# No.    Name      Ver    Type      Cards   Dimensions   Format
#   0  PRIMARY       1 PrimaryHDU       4   ()      
#   1  FIELDS        1 BinTableHDU     24   65534R x 5C   [1J, 1J, 1J, 16A, 1J]   
#   2  OBJECTS       1 BinTableHDU    186   39103516R x 59C   ['1J', '1E', '1E', '1E', '1E', '1J', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1E', '1E', '16A', '16A', '1I', '1E', '1E', '1E', '1J', '1I', '1I', '1E', '1E']   
# saved as /disks/shear10/ssli/K1000CS/skills_v07D7_Inputs/skills_v07D7_LF_321_kidsPhotometry_shear_noSG_newCut_A1_WeiCut.cat
# Elapsed:13:34.61,User=579.796,System=1470.965,CPU=251.7%.