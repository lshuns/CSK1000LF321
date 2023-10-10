# -*- coding: utf-8 -*-
# @Author: ssli
# @Date:   2022-11-13 15:43:33
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-22 17:16:31

### split catalogue into N and S patches
###### apply SOM gold selection
###### rename for easier use for cosmological analysis

import numpy as np
import pandas as pd 

from astropy.io import fits
from astropy.table import Table

# >>>>>>>>>>>>>>>>>>>>>>> general info
# Individual PATCH Designations
PATCHLIST = ["N", "S"]
# PATCH Lower RA limits
PATCHLO = [120, 315]
# PATCH Upper RA limits
PATCHHI = [240, 60]

# renaming
renames = {'RAJ2000': 'ALPHA_J2000', 'DECJ2000': 'DELTA_J2000',
            'AlphaRecalC_variance': 'AlphaRecal_variance', 
            'AlphaRecalC_weight': 'AlphaRecal_weight',
            'Flag_SOM_Fid_NONE': 'Flag_SOM_Fid',
            'AlphaRecalD2_e1': 'AlphaRecal_e1',
            'AlphaRecalD2_e2': 'AlphaRecal_e2'}

# >>>>>>>>>>>>>>>>>>>>>>> I/O

# the shear catalogue
inpath = '/disks/shear10/ssli/K1000CS/LF321_Outputs//K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A1_goldclasses.cat.A2.feather'

# where to save
outpath = '/disks/shear10/ssli/K1000CS/LF321_Outputs/K1000_{:}_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses.cat'

# the LDAC template
ldac_template = "../../utils/LDAC.cat"

# csv file includes wanted columns
wanted_cols = '../utils/wanted_cols_LF321_A12_SOM.csv'

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
## all float32 to float64
cata[cata.select_dtypes(np.float32).columns] = cata.select_dtypes(np.float32).astype(np.float64)
## rename columns
cata.rename(columns=renames, inplace=True)
## select gold sample
cata = cata[cata['Flag_SOM_Fid']>0]
cata.reset_index(drop=True, inplace=True)
print('number gold', len(cata))

# 3. loop over and cut
Nsum = 0
for i_p, patch in enumerate(PATCHLIST):

    TMPLO = PATCHLO[i_p]
    TMPHI = PATCHHI[i_p]

    # Normal RA limits: pick lo < RA <= hi
    if TMPLO < TMPHI:
        cata_selec = cata[(cata['ALPHA_J2000']>TMPLO)&(cata['ALPHA_J2000']<=TMPHI)]
    else:
        cata_selec = cata[(cata['ALPHA_J2000']>TMPLO)|(cata['ALPHA_J2000']<=TMPHI)]

    Ntmp = len(cata_selec)
    print('>>> number', patch, Ntmp)
    Nsum += Ntmp

    # to hdu
    ## build the col in hdul format
    cols_list = []
    for i_col, col in enumerate(col_names):
        cols_list.append(fits.Column(name=col, array=cata_selec[col].values, format=col_format[i_col]))
    del cata_selec
    hdu = fits.BinTableHDU.from_columns(cols_list)
    del cols_list
    hdu.name = "OBJECTS"
    ## add dummy TCOMM values
    for i, col in enumerate(col_names):
        hdu.header.set('TCOMM%d' % (i + 1), col,
                               after='TTYPE%d' % (i + 1))

    # combine with template LDAC HDU
    outpath_tmp = outpath.format(patch)
    with fits.open(ldac_template) as ldac:
        ldac.append(hdu)
        del hdu
        print("LDAC info:")
        ldac.info()
        ldac.writeto(outpath_tmp, overwrite=True)
    print(f"saved as {outpath_tmp}")

print('number in all patches', Nsum)

# (py377) veersemeer [129] > python H_split_NS.py 
# number original 23401764
# number gold 23401764
# >>> number N 10940568
# LDAC info:
# Filename: ../../utils/LDAC.cat
# No.    Name      Ver    Type      Cards   Dimensions   Format
#   0  PRIMARY       1 PrimaryHDU       4   ()      
#   1  FIELDS        1 BinTableHDU     24   65534R x 5C   [1J, 1J, 1J, 16A, 1J]   
#   2  OBJECTS       1 BinTableHDU    357   10940568R x 116C   ['16A', '1J', '1J', '1E', '1E', '1D', '1D', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1J', '1E', '1I', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1I', '1E', '1E']   
# saved as /disks/shear10/ssli/K1000CS/LF321_Outputs/K1000_N_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses.cat
# >>> number S 12461196
# LDAC info:
# Filename: ../../utils/LDAC.cat
# No.    Name      Ver    Type      Cards   Dimensions   Format
#   0  PRIMARY       1 PrimaryHDU       4   ()      
#   1  FIELDS        1 BinTableHDU     24   65534R x 5C   [1J, 1J, 1J, 16A, 1J]   
#   2  OBJECTS       1 BinTableHDU    357   12461196R x 116C   ['16A', '1J', '1J', '1E', '1E', '1D', '1D', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1J', '1E', '1I', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1E', '1E', '1E', '1E', '1I', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1I', '1E', '1E']   
# saved as /disks/shear10/ssli/K1000CS/LF321_Outputs/K1000_S_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses.cat
# number in all patches 23401764
# Elapsed:6:12.82,User=327.023,System=221.311,CPU=147.0%.