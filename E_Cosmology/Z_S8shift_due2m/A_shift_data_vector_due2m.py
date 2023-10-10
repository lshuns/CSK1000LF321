# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-05-06 11:28:08
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-18 22:27:43

### build the data vector with different m values

import os
import io
import sys

import numpy as np
import pandas as pd

from astropy.io import fits

sys.path.append('../external')
import wrapper_twopoint2 as wtp2

# >>>>>>>>>>>>>>>>> I/O

## fiducial m file
inpath_m = '../../B_Shear_Bias/results/m_skills_v07D7p1_LF_321_kidsPhotometry_A12_gold_Rewei_varCorr_dz0p1_2D_PSFmodellingCorr.csv'

## dm files
# dm_label_names = ['q', 'n']
dm_label_names = ['size']
variation_names = ['U', 'D']
indir_dm = '/net/grecht/data2/ssli_files/Projects/Projects/8ShearBias_ImSim/MultiBand_ImSim/sensitivity_test/galaxy/results'
dm_file_template = 'dm_ZBbins_skills_v07D7_{:}{:}_nogold_reweighted.csv'

## the theta options
theta_cosebis_list = ['theta_0.5_300', 'theta_2.0_300']
theta_cov_list = ['0.50_300.00', '2.00_300.00']

## the input info
input_cosebis_template = '/disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/COSEBIS/nmodes5/\
En_COSEBIS_K1000_ALL_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_{:}_zbins_{:}_{:}.asc'
N_Zbins = 5

## the covariance
covName_template = '/disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/COV/output/Covariance_blindNONE_nMaximum_5_{:}_nBins5.ascii'

## where to find NofZ
inpath_prefix = '/disks/shear10/ssli/K1000CS/LF321_Outputs/K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A1_DIRcols_Fid_blindNONE_TOMO'
nOfZNameList = [f'{inpath_prefix}{i_bin+1}_Nz.asc' for i_bin in range(5)]

## other statistics
### number density of galaxies per arcmin^2
filename = '/disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/COV/input/neff.txt'
nGalList = np.loadtxt(filename).tolist()
### sigma e
filename = '/disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/COV/input/sigmae.txt'
sigmaEpsList = np.loadtxt(filename).tolist()

## where to save
saveName_template  = '/disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/data_vector/cosebis_K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_goldclasses_Flag_SOM_Fid_filt_{:}_nbins5_theta_{:}_cosmosis.fits'

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

# load fiducial m 
cata_tmp = pd.read_csv(inpath_m)
## the first is for whole
m = (cata_tmp.loc[1:, 'm1'].values + cata_tmp.loc[1:, 'm2'].values)/2.
print('>>> m (fiducial)', m)
del cata_tmp

# loop over theta
for i_theta, theta_cosebis in enumerate(theta_cosebis_list):
    theta_cov = theta_cov_list[i_theta]
    covName = covName_template.format(theta_cov)
    theta_min, theta_max = map(float, theta_cov.split("_"))

    # loop over different dm
    for label in dm_label_names:
        for variation in variation_names:
            inpath_dm = os.path.join(indir_dm, dm_file_template.format(label, variation))
            cata_tmp = pd.read_csv(inpath_dm)
            ## the first is for whole
            dm = (cata_tmp.loc[1:, 'm1'].values + cata_tmp.loc[1:, 'm2'].values)/2.
            del cata_tmp
            m_test = m + dm

            print('=== ', label, variation)
            print('>>> dm', dm)
            print('>>> new m', m_test)

            # build m-corrected COSEBIs
            input_files = []
            m_corr_all  = []
            for bin1 in range(N_Zbins):
                for bin2 in range(bin1, N_Zbins):
                    fileNameInput = input_cosebis_template.format(theta_cosebis, bin1+1, bin2+1) 
                    input_files.append(fileNameInput)
                    m_corr = (1.+m_test[bin2])*(1.+m_test[bin1])
                    m_corr_all.append(m_corr)
            _, COSEBIs_vector_with_m_bias  = make_2pt_vector(input_files, np.asarray(m_corr_all))
            ## save to a tmp file in current dir
            tmp_file = './tmp_file_for_A_shift_data_vector_due2m.txt'
            np.savetxt(tmp_file, COSEBIs_vector_with_m_bias)

            # combine to a data vector
            saveName = saveName_template.format(f'test_m_{label}{variation}', theta_cov)
            wtp2.saveFitsTwoPoint(
                    nbTomoN=0, nbTomoG=5,
                    N_theta=9, theta_min=theta_min, theta_max=theta_max,
                    N_ell=8, ell_min=100, ell_max=1500,
                    nbModes=5,
                    prefix_Flinc=None,
                    prefix_CosmoSIS=None,
                    scDict={'use_stats': 'En'.lower()},
                    meanTag='file', meanName=tmp_file,
                    covTag='file', covName=covName,
                    nOfZNameList=nOfZNameList, nGalList=nGalList, sigmaEpsList=sigmaEpsList,
                    saveName=saveName)
            os.remove(tmp_file)

# >>> m (fiducial) [-0.02067462 -0.02324236 -0.0153159   0.01488055  0.03133205]
# ===  q U
# >>> dm [-0.0003   0.0022   0.00125  0.0036   0.00635]
# >>> new m [-0.02097462 -0.02104236 -0.0140659   0.01848055  0.03768205]
# Saved "/disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/data_vector/cosebis_K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_goldclasses_Flag_SOM_Fid_filt_test_m_qU_nbins5_theta_0.50_300.00_cosmosis.fits"
# ===  q D
# >>> dm [-0.0016  -0.0027  -0.00105 -0.00405 -0.0028 ]
# >>> new m [-0.02227462 -0.02594236 -0.0163659   0.01083055  0.02853205]
# Saved "/disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/data_vector/cosebis_K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_goldclasses_Flag_SOM_Fid_filt_test_m_qD_nbins5_theta_0.50_300.00_cosmosis.fits"
# ===  n U
# >>> dm [-0.00155 -0.0027  -0.0012  -0.00165 -0.0017 ]
# >>> new m [-0.02222462 -0.02594236 -0.0165159   0.01323055  0.02963205]
# Saved "/disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/data_vector/cosebis_K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_goldclasses_Flag_SOM_Fid_filt_test_m_nU_nbins5_theta_0.50_300.00_cosmosis.fits"
# ===  n D
# >>> dm [0.00065 0.0006  0.0022  0.0017  0.00315]
# >>> new m [-0.02002462 -0.02264236 -0.0131159   0.01658055  0.03448205]
# Saved "/disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/data_vector/cosebis_K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_goldclasses_Flag_SOM_Fid_filt_test_m_nD_nbins5_theta_0.50_300.00_cosmosis.fits"
# ===  q U
# >>> dm [-0.0003   0.0022   0.00125  0.0036   0.00635]
# >>> new m [-0.02097462 -0.02104236 -0.0140659   0.01848055  0.03768205]
# Saved "/disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/data_vector/cosebis_K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_goldclasses_Flag_SOM_Fid_filt_test_m_qU_nbins5_theta_2.00_300.00_cosmosis.fits"
# ===  q D
# >>> dm [-0.0016  -0.0027  -0.00105 -0.00405 -0.0028 ]
# >>> new m [-0.02227462 -0.02594236 -0.0163659   0.01083055  0.02853205]
# Saved "/disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/data_vector/cosebis_K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_goldclasses_Flag_SOM_Fid_filt_test_m_qD_nbins5_theta_2.00_300.00_cosmosis.fits"
# ===  n U
# >>> dm [-0.00155 -0.0027  -0.0012  -0.00165 -0.0017 ]
# >>> new m [-0.02222462 -0.02594236 -0.0165159   0.01323055  0.02963205]
# Saved "/disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/data_vector/cosebis_K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_goldclasses_Flag_SOM_Fid_filt_test_m_nU_nbins5_theta_2.00_300.00_cosmosis.fits"
# ===  n D
# >>> dm [0.00065 0.0006  0.0022  0.0017  0.00315]
# >>> new m [-0.02002462 -0.02264236 -0.0131159   0.01658055  0.03448205]
# Saved "/disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/data_vector/cosebis_K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_goldclasses_Flag_SOM_Fid_filt_test_m_nD_nbins5_theta_2.00_300.00_cosmosis.fits"
