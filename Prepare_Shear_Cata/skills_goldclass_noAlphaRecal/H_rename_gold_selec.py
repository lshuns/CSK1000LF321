# -*- coding: utf-8 -*-
# @Author: ssli
# @Date:   2022-11-13 15:43:33
# @Last Modified by:   ssli
# @Last Modified time: 2022-12-25 12:23:13

### apply SOM gold selection 
###### and rename for consistent with KiDS catalogue

import numpy as np
import pandas as pd 

from astropy.io import fits
from astropy.table import Table

# >>>>>>>>>>>>>>>>>>>>>>> general info

# renaming
renames = {'X_WORLD': 'ALPHA_J2000', 'Y_WORLD': 'DELTA_J2000',
            'e1_LF_r': 'autocal_e1', 'e2_LF_r': 'autocal_e2',
            'SNR_LF_r': 'model_SNratio', 
            'scalelength_LF_r': 'autocal_scalelength_pixels',
            'oldweight_LF_r': 'weight',
            'psf_e1_LF_r': 'PSF_e1', 'psf_e2_LF_r': 'PSF_e2',
            'psf_strehl_LF_r': 'PSF_Strehl_ratio', 
            'psf_Q11_LF_r': 'PSF_Q11', 'psf_Q22_LF_r': 'PSF_Q22', 'psf_Q12_LF_r': 'PSF_Q12',
            'class_LF_r': 'fitclass', 'contamination_radius_LF_r': 'contamination_radius',
            'nm_LF_r': 'neighbour_mag', 'nd_LF_r': 'neighbour_distance',
            'LS_variance_LF_r': '2D_measurement_variance',
            'Flag_SOM_Fid_NONE': 'Flag_SOM_Fid'}

# >>>>>>>>>>>>>>>>>>>>>>> I/O

# the shear catalogue
inpath = '/disks/shear10/ssli/K1000CS//skills_v07D7_noA_Inputs//skills_v07D7_LF_321_kidsPhotometry_shear_noSG_newCut_WeiCut_goldclasses.cat'

# where to save
outpath = inpath.replace('.cat', '_goldSelected.cat')

# the LDAC template
ldac_template = "../../utils/LDAC.cat"

# csv file includes wanted columns
wanted_cols = '../utils/wanted_cols_skills_v07D7_noA_SOM.csv'

# >>>>>>>>>>>>>>>>>>>>>>>> workhorse

# 1. column info
tmp = pd.read_csv(wanted_cols)
col_names = tmp['col'].to_list()
col_format = tmp['format'].to_list()
del tmp
## rename
for key, value in renames.items():
    for i_col, col_name in enumerate(col_names):
        if col_name == key:
            col_names[i_col] = value

# 2. get catalogues
file_type = inpath[-3:]
if file_type == 'her':
    cata = pd.read_feather(inpath)
elif file_type == 'cat':
    with fits.open(inpath) as hdul:
        cata = Table(hdul['OBJECTS'].data).to_pandas()
elif file_type == 'its':
    with fits.open(inpath) as hdul:
        cata = Table(hdul[1].data).to_pandas()
else:
    raise Exception(f'Not supported input file type! {inpath}')
print('number original', len(cata))
## rename columns
cata.rename(columns=renames, inplace=True)
## select gold sample
cata = cata[cata['Flag_SOM_Fid']>0]
cata.reset_index(drop=True, inplace=True)
print('number gold', len(cata))

# 3. to LDAC cat
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

# combine with template LDAC HDU
with fits.open(ldac_template) as ldac:
    ldac.append(hdu)
    del hdu
    print("LDAC info:")
    ldac.info()
    ldac.writeto(outpath, overwrite=True)
print(f"saved as {outpath}")

# number original 39,251,930
# number gold 29,253,262
# LDAC info:
# Filename: ../../utils/LDAC.cat
# No.    Name      Ver    Type      Cards   Dimensions   Format
#   0  PRIMARY       1 PrimaryHDU       4   ()      
#   1  FIELDS        1 BinTableHDU     24   65534R x 5C   [1J, 1J, 1J, 16A, 1J]   
#   2  OBJECTS       1 BinTableHDU    174   29253262R x 55C   ['1J', '1E', '1E', '1E', '1E', '1J', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1E', '1E', '16A', '16A', '1I', '1E', '1E', '1E', '1I']   
# saved as /disks/shear10/ssli/K1000CS//skills_v07D7_noA_Inputs//skills_v07D7_LF_321_kidsPhotometry_shear_noSG_newCut_WeiCut_goldclasses_goldSelected.cat
# Elapsed:9:51.29,User=365.787,System=161.019,CPU=89.0%.
