# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-04-15 16:17:33
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-06 12:33:07

### combine the individual COSEBIs 
###### This is a hack from https://github.com/KiDS-WL/Cat_to_Obs_K1000_P1/blob/master/2pt_data_to_fits/MakeDataVectors.py
######### We only preserve those related to the COSEBIS

import os

import numpy as np
import pandas as pd

# >>>>>>>>>>>>>>>>>>> I/O

# the input info
# theta_tag = 'theta_0.5_300'
theta_tag = 'theta_2.0_300'
input_prefix = f'/disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/COSEBIS/nmodes5/\
En_COSEBIS_K1000_ALL_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_{theta_tag}_zbins'
N_Zbins = 5

# where to save
outdir = '/disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/COSEBIS/nmodes5/'

# the fiducial shear bias
inpath_m = '../../B_Shear_Bias/results/m_skills_v07D7p1_LF_321_kidsPhotometry_A12_gold_Rewei_varCorr_dz0p1_2D_PSFmodellingCorr.csv'
cata_tmp = pd.read_csv(inpath_m)
## the first is for whole
m = (cata_tmp.loc[1:, 'm1'].values + cata_tmp.loc[1:, 'm2'].values)/2.
print('>>> m (fiducial)', m)
del cata_tmp

# the raw shear bias
inpath_m = '../../B_Shear_Bias/results/m_skills_v07D7p1_LF_321_kidsPhotometry_A12_gold_Rewei.csv'
cata_tmp = pd.read_csv(inpath_m)
## the first is for whole
m_raw = (cata_tmp.loc[1:, 'm1'].values + cata_tmp.loc[1:, 'm2'].values)/2.
print('>>> m (raw)', m_raw)
del cata_tmp

# ## outinfo
# >>> m (fiducial) [-0.02067462 -0.02324236 -0.0153159   0.01488055  0.03133205]
# >>> m (raw) [-0.02335 -0.02485 -0.0132   0.01795  0.0323 ]
# >>> raw COSEBIs saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/COSEBIS/nmodes5/cosebis_combined_theta_2.0_300_no_m_bias.asc
# >>> fiducial m corrected COSEBIs saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/COSEBIS/nmodes5/cosebis_combined_theta_2.0_300_with_m_bias.asc
# >>> raw m corrected COSEBIs saved to /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/COSEBIS/nmodes5/cosebis_combined_theta_2.0_300_with_raw_m_bias.asc

# >>>>>>>>>>>>>>>>>>> the functions
# Reads in from the list of input_files and puts them all into a long vector. 
# Make sure that the ordering is correct, col starts from 1 instead of 0
def make_2pt_vector(input_files, m_corr, col=1):
    for rp in range(len(input_files)):
        file= open(input_files[rp])
        data=np.loadtxt(file,comments='#')
        if data.ndim==1:
            if rp==0:
                data_all      = data.copy()
                data_all_corr = data/m_corr[rp]
            else:
                data_all      = np.hstack((data_all,data))
                data_all_corr = np.hstack((data_all_corr,data/m_corr[rp]))
        else:
            if rp==0:
                data_all      = data[:,col-1].copy()
                data_all_corr = data[:,col-1]/m_corr[rp]
            else:
                data_all      = np.hstack((data_all,data[:,col-1]))
                data_all_corr = np.hstack((data_all_corr,data[:,col-1]/m_corr[rp]))
    return data_all,data_all_corr

# >>>>>>>>>>>>>>>>>> the workhorse

# COSEBIs
input_files = []
m_corr_all  = []
m_corr_all_raw = []
for bin1 in range(N_Zbins):
    for bin2 in range(bin1, N_Zbins):
        ## the measurements
        fileNameInput = f'{input_prefix}_{bin1+1}_{bin2+1}.asc'
        input_files.append(fileNameInput)
        ## fiducial m
        m_corr_all.append((1.+m[bin2])*(1.+m[bin1]))
        ## raw m
        m_corr_all_raw.append((1.+m_raw[bin2])*(1.+m_raw[bin1]))

COSEBIs_vector_no_m_bias, COSEBIs_vector_with_m_bias  = make_2pt_vector(input_files, np.asarray(m_corr_all))
_, COSEBIs_vector_with_raw_m_bias  = make_2pt_vector(input_files, np.asarray(m_corr_all_raw))

name_tag = theta_tag + '_no_m_bias'
savename = os.path.join(outdir, f'cosebis_combined_{name_tag}.asc')
np.savetxt(savename,COSEBIs_vector_no_m_bias)
print('>>> raw COSEBIs saved to', savename)

name_tag = theta_tag + '_with_m_bias'
savename = os.path.join(outdir, f'cosebis_combined_{name_tag}.asc')
np.savetxt(savename,COSEBIs_vector_with_m_bias)
print('>>> fiducial m corrected COSEBIs saved to', savename)

name_tag = theta_tag + '_with_raw_m_bias'
savename = os.path.join(outdir, f'cosebis_combined_{name_tag}.asc')
np.savetxt(savename,COSEBIs_vector_with_raw_m_bias)
print('>>> raw m corrected COSEBIs saved to', savename)