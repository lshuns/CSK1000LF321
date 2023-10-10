# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-04-12 17:41:36
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-22 17:23:14

### cut catalogues into tomographic bins
#### and apply c-term correction
###### This is a hack from https://github.com/KiDS-WL/Cat_to_Obs_K1000_P1/blob/master/Calc_2pt_Stats/doall_calc2pt.sh
######### mode: \"CREATETOMO\": cut catalogues into tomographic bins and calculate and subtract c-term"

import os
import sys
import ldac
import numpy as np
from astropy.io import fits

# >>>>>>>>> I/O

# the input info
indir = '/disks/shear10/ssli/K1000CS/LF321_Outputs/'
file_tag = 'V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses'
patches = ['K1000_N', 'K1000_S']

# the out info
outdir = '/disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/TOMOCATS/'
save_tag = 'BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid'

# column info
e1_col = 'AlphaRecal_e1'
e2_col = 'AlphaRecal_e2'
wt_col = 'AlphaRecal_weight'
zb_col = 'Z_B'
psfe1_col = 'PSF_e1'
psfe2_col = 'PSF_e2'
ra_col = 'ALPHA_J2000'
dec_col = 'DELTA_J2000'
x_col = 'Xpos'
y_col = 'Ypos'
gold_col = 'Flag_SOM_Fid'

# the delta_epsf(X,Y) data
## Note we used to call this "c" - but it's not a c-term
## It's a delta_epsf term
c1_map_file = '../PSFRES_CORRMAP/c1_map.fits'
c2_map_file = '../PSFRES_CORRMAP/c2_map.fits'
exp_map_file = '../PSFRES_CORRMAP/exposure_map.fits'

# the binning range
N_Zbins = 5
ZB_ranges = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2]

# >>>>>>>>> workhorse

# Read in the delta_epsf(X,Y) data
with fits.open(c1_map_file) as hdul:
    depsf_1_map = hdul[0].data
with fits.open(c2_map_file) as hdul:
    depsf_2_map = hdul[0].data

# find the edges of the image from the exposure map
with fits.open(exp_map_file) as hdul:
    exp_map = hdul[0].data
# and clip so it is zero for out of the image and 1 within
w_map = np.clip(exp_map,0,1)
w_map = w_map.astype(int)

# remove average depsf1/depsf2
ave1 = np.sum(w_map*depsf_1_map)/np.sum(w_map)
ave2 = np.sum(w_map*depsf_2_map)/np.sum(w_map)

depsf_1_map = (depsf_1_map - ave1)*w_map
depsf_2_map = (depsf_2_map - ave2)*w_map

# loop over all inputs
for patch in patches:

    # open the ldac catalogue using functions in ldac.py
    infile = os.path.join(indir, f'{patch}_{file_tag}.cat')
    ldac_cat = ldac.LDACCat(infile)
    ldac_table = ldac_cat['OBJECTS']
    del ldac_cat

    # needed values
    e1 = ldac_table[e1_col]
    e2 = ldac_table[e2_col]
    weight = ldac_table[wt_col]
    Z_B = ldac_table[zb_col]
    PSF_e1 = ldac_table[psfe1_col]
    PSF_e2 = ldac_table[psfe2_col]
    ALPHA_J2000 = ldac_table[ra_col]
    DELTA_J2000 = ldac_table[dec_col]
    Xpos_in = ldac_table[x_col]
    Ypos_in = ldac_table[y_col]
    FLAG_SOM = ldac_table[gold_col]
    del ldac_table

    # calculate the delta_epsf(x,y) term
    Nobj = np.shape(Xpos_in)[0]
    XYpos = np.vstack((Xpos_in,Ypos_in))
    XYpos_short = XYpos/10.
    del XYpos
    XYpos_int = XYpos_short.astype(int)
    # mock delta_epsf2 and delta_epsf2 columns
    depsf_xy_1 = np.zeros(Nobj)
    depsf_xy_2 = np.zeros(Nobj)
    for i in range(Nobj):
        depsf_xy_1[i] = depsf_1_map[XYpos_int[0,i],XYpos_int[1,i]]
        depsf_xy_2[i] = depsf_2_map[XYpos_int[0,i],XYpos_int[1,i]]
    del XYpos_int

    # split in tomo bins
    for i_zbin in range(N_Zbins):

        outfile = os.path.join(outdir, f'{patch}_{save_tag}_{N_Zbins}Z_{i_zbin+1}.fits')
        zmin = ZB_ranges[i_zbin]
        zmax = ZB_ranges[i_zbin+1]
        print(f'>>> selected Z range ({zmin}, {zmax}]')

        ztomo = ((Z_B<=zmax) & (Z_B>zmin) & (FLAG_SOM > 0))
        # apply the tomographic/SOM selection
        e1_inbin = e1[ztomo]
        e2_inbin = e2[ztomo]
        ra_inbin = ALPHA_J2000[ztomo]
        dec_inbin = DELTA_J2000[ztomo]
        PSF_e1_inbin = PSF_e1[ztomo]
        PSF_e2_inbin = PSF_e2[ztomo]
        w_inbin = weight[ztomo]
        depsf_xy_1_inbin = depsf_xy_1[ztomo]
        depsf_xy_2_inbin = depsf_xy_2[ztomo]
        del ztomo

        # carry through the square of the weight for
        # Npair calculation hack with Treecorr
        wsq_inbin = w_inbin*w_inbin

        # weighted mean   
        c1 = np.average(e1_inbin, weights=w_inbin)
        c2 = np.average(e2_inbin, weights=w_inbin)
        # apply correction
        e1_corr = e1_inbin - c1
        e2_corr = e2_inbin - c2
        del c1, c2, e1_inbin, e2_inbin

        # Write out to output file - crucial that RA/DEC (in degrees) are double precision
        # If you don't have that you round to a couple of arcsec for fields with ra > 100
        hdulist = fits.BinTableHDU.from_columns(
            [fits.Column(name='ALPHA_J2000', format='1D', unit='deg',array=ra_inbin),
             fits.Column(name='DELTA_J2000', format='1D', unit='deg',array=dec_inbin),
             fits.Column(name='e1', format='1E', array=e1_corr),
             fits.Column(name='e2', format='1E', array=e2_corr),
             fits.Column(name='PSF_e1', format='1E', array=PSF_e1_inbin),
             fits.Column(name='PSF_e2', format='1E', array=PSF_e2_inbin),
             fits.Column(name='dPSF_e1_xy', format='1E', array=depsf_xy_1_inbin),
             fits.Column(name='dPSF_e2_xy', format='1E', array=depsf_xy_2_inbin),
             fits.Column(name='weight', format='1E', array=w_inbin),
             fits.Column(name='weightsq', format='1E', array=wsq_inbin)])
        del ra_inbin, dec_inbin, e1_corr, e2_corr, PSF_e1_inbin, PSF_e2_inbin
        del depsf_xy_1_inbin, depsf_xy_2_inbin, w_inbin, wsq_inbin 
        hdulist.writeto(outfile, overwrite=True)
        del hdulist
        print('saved to', outfile)

# (py377) veersemeer [133] > python A_CREATETOMO.py 
# >>> selected Z range (0.1, 0.3]
# saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/TOMOCATS/K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_5Z_1.fits
# >>> selected Z range (0.3, 0.5]
# saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/TOMOCATS/K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_5Z_2.fits
# >>> selected Z range (0.5, 0.7]
# saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/TOMOCATS/K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_5Z_3.fits
# >>> selected Z range (0.7, 0.9]
# saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/TOMOCATS/K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_5Z_4.fits
# >>> selected Z range (0.9, 1.2]
# saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/TOMOCATS/K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_5Z_5.fits
# >>> selected Z range (0.1, 0.3]
# saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/TOMOCATS/K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_5Z_1.fits
# >>> selected Z range (0.3, 0.5]
# saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/TOMOCATS/K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_5Z_2.fits
# >>> selected Z range (0.5, 0.7]
# saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/TOMOCATS/K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_5Z_3.fits
# >>> selected Z range (0.7, 0.9]
# saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/TOMOCATS/K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_5Z_4.fits
# >>> selected Z range (0.9, 1.2]
# saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/TOMOCATS/K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_5Z_5.fits
# Elapsed:1:50.56,User=61.380,System=11.067,CPU=65.5%.
