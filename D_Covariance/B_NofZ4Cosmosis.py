# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-04-18 18:03:05
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-24 10:42:33

### prepare NofZ for cosmoSIS
###### a hack of https://github.com/AngusWright/CosmoPipe/blob/master/scripts/MakeNofZForCosmosis_function.py

import os

import numpy as np

from astropy.io import fits

# >>>>>>>>>>>>>>>>> I/O

# where to find the separate NofZ files
indir = '/disks/shear10/ssli/K1000CS/LF321_Outputs'
file_prefix = 'K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A1_DIRcols_Fid_blindNONE_TOMO'
file_suffix = '_Nz.asc' 

# where to save
outdir = '/disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/COV/input'

# number of tomo bins
nBins = 5

# >>>>>>>>>>>>>>>>> workhorse

# get the NofZ files
input_files = [os.path.join(indir, file_prefix + str(bin1+1) + file_suffix) for bin1 in range(nBins)]
print('number of tomo bins', nBins)

# combine values
cols = []
for bin1, input_file in enumerate(input_files):

    nofZ = np.loadtxt(input_file)

    if(bin1==0):
        DeltaZ = nofZ[1, 0] - nofZ[0, 0]
        Z_LOW = nofZ[:, 0]
        Z_HIGH = nofZ[:, 0] + DeltaZ
        Z_MID = Z_LOW + DeltaZ/2.

        cols.append(fits.Column(name='Z_lOW', format='D', array=Z_LOW))
        cols.append(fits.Column(name='Z_HIGH', format='D', array=Z_HIGH))
        cols.append(fits.Column(name='Z_MID', format='D', array=Z_MID))

    cols.append(fits.Column(name='BIN'+str(bin1+1), format='D', array=nofZ[:, 1]))

# build a hdul and save
outpath = os.path.join(outdir, file_prefix + 'comb_Nz.fits')
new_cols = fits.ColDefs(cols)
hdulist_new = fits.BinTableHDU.from_columns(new_cols)
hdulist_new.header['NZDATA'] = True
hdulist_new.header['EXTNAME'] = 'NZ_source'
hdulist_new.header['NBIN'] = nBins
hdulist_new.header['NZ'] = len(Z_LOW)
hdulist_new.writeto(outpath)
print('saved to', outpath)
