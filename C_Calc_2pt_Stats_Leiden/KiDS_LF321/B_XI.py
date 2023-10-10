# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-04-12 19:32:12
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-15 11:58:58

### calculate the two-point correlation functions
###### This is a hack from https://github.com/KiDS-WL/Cat_to_Obs_K1000_P1/blob/master/Calc_2pt_Stats/doall_calc2pt.sh
######### mode: \"XI\": calculate xi_+/- for tomo bin pair i j"

import os
import sys
import treecorr
import numpy as np
from astropy.io import fits

# >>>>>>>>>>>>>>>> I/O

# general setups
# nbins = 4000 
nbins = 9
theta_min = 0.5 # arcmin
# theta_min = 2 # arcmin
theta_max = 300.0 # arcmin
num_threads = 20 # number of cores

# where to save
outdir = '/disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI'

# where to find the inputs
indir = '/disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/TOMOCATS'
file_tag = 'BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid'
patches = ['K1000_N', 'K1000_S']
N_Zbins = 5

# >>>>>>>>>>>>>>>>> functions
######## gen_write from old treecorr
def gen_write(file_name, col_names, columns, params=None, precision=4, file_type=None, logger=None):
    """Write some columns to an output file with the given column names.
    We do this basic functionality a lot, so put the code to do it in one place.
    :param file_name:   The name of the file to write to.
    :param col_names:   A list of columns names for the given columns.
    :param columns:     A list of numpy arrays with the data to write.
    :param params:      A dict of extra parameters to write at the top of the output file (for
                        ASCII output) or in the header (for FITS output).  (default: None)
    :param precision:   Output precision for ASCII. (default: 4)
    :param file_type:   Which kind of file to write to. (default: determine from the file_name
                        extension)
    :param logger:      If desired, a logger object for logging. (default: None)
    """
    if len(col_names) != len(columns):
        raise ValueError("col_names and columns are not the same length.")
    if len(columns) == 0:
        raise ValueError("len(columns) == 0")
    for col in columns[1:]:
        if col.shape != columns[0].shape:
            raise ValueError("columns are not all the same shape")

    columns = [ col.flatten() for col in columns ]
    writer = treecorr.util.make_writer(file_name, precision, file_type, logger)
    with writer:
        writer.write(col_names, columns, params=params)

# >>>>>>>>>>>>>>>>> workhorse

# calculate XI for each patch
for patch in patches:
    print('>>> working on', patch)

    for i_bin1 in range(N_Zbins):
        for i_bin2 in range(i_bin1, N_Zbins):

            # the files
            outfile = os.path.join(outdir, f'XI_{patch}_{file_tag}_nbins_{nbins}_theta_{theta_min:.1f}_{theta_max:.1f}_zbins_{i_bin1+1}_{i_bin2+1}.asc')
            fitscat1 = os.path.join(indir, f'{patch}_{file_tag}_{N_Zbins}Z_{i_bin1+1}.fits')
            fitscat2 = os.path.join(indir, f'{patch}_{file_tag}_{N_Zbins}Z_{i_bin2+1}.fits')

            # construct treecorr catalogues
            cat1 = treecorr.Catalog(fitscat1, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                              g1_col='e1', g2_col='e2', w_col='weight')
            cat2 = treecorr.Catalog(fitscat2, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                              g1_col='e1', g2_col='e2', w_col='weight')

            if nbins > 100: ## Fine-binning
                inbinslop = 1.5
            else: ## Broad bins
                inbinslop = 0.08

            # Define the binning based on command line input
            # Log is the default bin_type for Treecorr
            gg = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
                        bin_slop=inbinslop)

            # Calculate the 2pt correlation function
            gg.process(cat1, cat2, num_threads=num_threads)
            del cat1, cat2

            # prepare the weighted_square catalogues - hack so that Treecorr returns the correct Npairs for a weighted sample
            cat1_wsq = treecorr.Catalog(fitscat1, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                          g1_col='e1', g2_col='e2', w_col='weightsq')
            cat2_wsq = treecorr.Catalog(fitscat2, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                          g1_col='e1', g2_col='e2', w_col='weightsq')

            # Define the binning based on command line input
            # Log is the default bin_type for Treecorr
            gg_wsq = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', \
                                            bin_slop=inbinslop)    

            # Calculate the weighted square 2pt correlation function
            gg_wsq.process(cat1_wsq, cat2_wsq, num_threads=num_threads)
            del cat1_wsq, cat2_wsq

            # Calculate the weighted Npairs = sum(weight_a*weight_b)^2 / sum(weight_a^2*weight_b^2)
            npairs_weighted = (gg.weight)*(gg.weight)/gg_wsq.weight

            #Use treecorr to write out the output file updating the npairs column and praise-be for Jarvis and his well documented code
            #as sigma_xip = sigma_xim, I've replaced sigma_xim with the raw npairs so we can store it in case useful at any point
            gen_write(outfile,
                    ['r_nom', 'meanr', 'meanlogr', 'xip', 'xim', 'xip_im', 'xim_im', 
                    'sigma_xip', 'npairs', 
                    'weight', 'npairs_weighted'],
                    [ gg.rnom, gg.meanr, gg.meanlogr, gg.xip, gg.xim, gg.xip_im, gg.xim_im, 
                    np.sqrt(gg.varxip), gg.npairs, 
                    gg.weight, npairs_weighted], 
                    precision=12)
            print('results saved to', outfile)

# >>> working on K1000_N
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_1_1.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_1_2.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_1_3.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_1_4.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_1_5.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_2_2.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_2_3.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_2_4.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_2_5.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_3_3.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_3_4.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_3_5.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_4_4.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_4_5.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_5_5.asc
# >>> working on K1000_S
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_1_1.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_1_2.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_1_3.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_1_4.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_1_5.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_2_2.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_2_3.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_2_4.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_2_5.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_3_3.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_3_4.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_3_5.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_4_4.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_4_5.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_5_5.asc
# Elapsed:46:40:37.87,User=3251670.237,System=5697.386,CPU=1938.4%.

# >>> working on K1000_N
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_1_1.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_1_2.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_1_3.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_1_4.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_1_5.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_2_2.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_2_3.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_2_4.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_2_5.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_3_3.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_3_4.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_3_5.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_4_4.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_4_5.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_5_5.asc
# >>> working on K1000_S
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_1_1.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_1_2.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_1_3.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_1_4.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_1_5.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_2_2.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_2_3.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_2_4.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_2_5.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_3_3.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_3_4.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_3_5.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_4_4.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_4_5.asc
# results saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins_5_5.asc
# Elapsed:19:55:40.04,User=4168381.244,System=7803.295,CPU=5821.2%.
