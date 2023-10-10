# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-04-25 09:27:30
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-06 12:04:19

### prepare the input fits for CosmoSIS
###### a hack of https://github.com/KiDS-WL/Cat_to_Obs_K1000_P1/blob/master/2pt_data_to_fits/save_and_check_Phase1.py

import os
import sys

import numpy as np

from astropy.io import fits

sys.path.append('./external')
import wrapper_twopoint2 as wtp2

# >>>>>>>>>>>>>>>>> I/O

## where to find the cosebis
cosebis_prefix = '/disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/COSEBIS/nmodes5/cosebis_combined'
cosebis_suffix = '.asc'
## where to find covariance
cov_prefix = '/disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/COV/output/Covariance_blindNONE_nMaximum_5'
cov_suffix = 'nBins5.ascii'
## where to save
saveName_template = '/disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/data_vector/cosebis_K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_goldclasses_Flag_SOM_Fid_filt_{:}_nbins5_theta_{:}_cosmosis.fits'

## the theta options
theta_cosebis_list = ['theta_0.5_300', 'theta_2.0_300']
theta_cov_list = ['0.50_300.00', '2.00_300.00']

## the m options
m_list = ['no_m_bias', 'with_raw_m_bias', 'with_m_bias']

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

# >>>>>>>>>>>>>>>>> workhorse

## save for COSEBIs
for i_theta, theta_cosebis in enumerate(theta_cosebis_list):
    theta_cov = theta_cov_list[i_theta]
    covName = f'{cov_prefix}_{theta_cov}_{cov_suffix}'
    theta_min, theta_max = map(float, theta_cov.split("_"))
    # print(">>> theta", theta_min, theta_max)
    for m_tag in m_list:
        cosebis_filename = f'{cosebis_prefix}_{theta_cosebis}_{m_tag}{cosebis_suffix}'
        saveName  = saveName_template.format(m_tag, theta_cov)

        wtp2.saveFitsTwoPoint(
                nbTomoN=0, nbTomoG=5,
                N_theta=9, theta_min=theta_min, theta_max=theta_max,
                N_ell=8, ell_min=100, ell_max=1500,
                nbModes=5,
                prefix_Flinc=None,
                prefix_CosmoSIS=None,
                scDict={'use_stats': 'En'.lower()},
                meanTag='file', meanName=cosebis_filename,
                covTag='file', covName=covName,
                nOfZNameList=nOfZNameList, nGalList=nGalList, sigmaEpsList=sigmaEpsList,
                saveName=saveName)