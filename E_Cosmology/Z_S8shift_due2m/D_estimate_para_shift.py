# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-05-22 11:37:17
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-22 12:54:07

### get the uncertainties introduced by the residual m

import os

import numpy as np
import pandas as pd

# ++++++++++++++++++++++ I/O
# run_name = 'multinest_nIA'
# err_range = [-0.023, 0.027]

run_name = 'polychord_nIA'
m_variations = ['test_m_sizeU', 'test_m_sizeD',
                'test_m_qU', 'test_m_qD',
                'test_m_nU', 'test_m_nD']
# ### out info

indir = f'/disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/{run_name}'
# results with the fiducial m
infile0 = f'maxpost_{run_name}_start.txt'
# results for comparison
infile_list = [f'maxpost_{run_name}_{m_variation}.txt' for m_variation in m_variations]

paras = ['s8', 'omega_m', 'sigma_8', 'a_ia']

# >>>>>>>>>>>>>>>>>>>>> the parameter names in MCMC chain
input_names={ 'omch2' :'cosmological_parameters--omch2', 
              'ombh2':'cosmological_parameters--ombh2', 
              'h'    :'cosmological_parameters--h0', 
              'h_out'     :'COSMOLOGICAL_PARAMETERS--H0',
              'n_s'   :'cosmological_parameters--n_s', 
              's8_in':'cosmological_parameters--s_8_input', 
              'logt_agn':'halo_model_parameters--logt_agn',
              'a_bar'      :'halo_model_parameters--a', 
              'a_ia':'intrinsic_alignment_parameters--a', 
              'alpha_ia': 'intrinsic_alignment_parameters--alpha',
              'deltaz_uncorr_1' :'nofz_shifts--uncorr_bias_1', 
              'deltaz_uncorr_2' :'nofz_shifts--uncorr_bias_2', 
              'deltaz_uncorr_3' :'nofz_shifts--uncorr_bias_3', 
              'deltaz_uncorr_4' :'nofz_shifts--uncorr_bias_4', 
              'deltaz_uncorr_5' :'nofz_shifts--uncorr_bias_5', 
              'uncorr_m1' :'shear_calibration_parameters--uncorr_m1', 
              'uncorr_m2' :'shear_calibration_parameters--uncorr_m2', 
              'uncorr_m3' :'shear_calibration_parameters--uncorr_m3', 
              'uncorr_m4' :'shear_calibration_parameters--uncorr_m4', 
              'uncorr_m5' :'shear_calibration_parameters--uncorr_m5', 
              'delta_c' : 'shear_c_bias--delta_c',
              's8':    'COSMOLOGICAL_PARAMETERS--S_8', 
              'sigma_8':  'COSMOLOGICAL_PARAMETERS--SIGMA_8', 
              'a_s': 'COSMOLOGICAL_PARAMETERS--A_S', 
              'omega_m':    'COSMOLOGICAL_PARAMETERS--OMEGA_M', 
              'omega_nu':    'COSMOLOGICAL_PARAMETERS--OMEGA_NU', 
              'omega_lambda':       'COSMOLOGICAL_PARAMETERS--OMEGA_LAMBDA', 
              'theta_mc':      'COSMOLOGICAL_PARAMETERS--COSMOMC_THETA', 
              'deltaz_1':'NOFZ_SHIFTS--BIAS_1', 
              'deltaz_2':'NOFZ_SHIFTS--BIAS_2', 
              'deltaz_3':'NOFZ_SHIFTS--BIAS_3', 
              'deltaz_4':'NOFZ_SHIFTS--BIAS_4', 
              'deltaz_5':'NOFZ_SHIFTS--BIAS_5', 
              'deltaz_out_1':'DELTA_Z_OUT--BIN_1', 
              'deltaz_out_2':'DELTA_Z_OUT--BIN_2', 
              'deltaz_out_3':'DELTA_Z_OUT--BIN_3', 
              'deltaz_out_4':'DELTA_Z_OUT--BIN_4', 
              'deltaz_out_5':'DELTA_Z_OUT--BIN_5',
              'm1':'SHEAR_CALIBRATION_PARAMETERS--M1',
              'm2':'SHEAR_CALIBRATION_PARAMETERS--M2',
              'm3':'SHEAR_CALIBRATION_PARAMETERS--M3',
              'm4':'SHEAR_CALIBRATION_PARAMETERS--M4',
              'm5':'SHEAR_CALIBRATION_PARAMETERS--M5',
              'prior':'prior', 
              'like':'like', 
              'post':'post', 
              'weight':'weight'}

# >>>>>>>>>>>>>>>>>>>>> workhorse

# load the fiducial results
inpath = os.path.join(indir, infile0)
## get column names
with open(inpath) as tmp:
    cols_names = (tmp.readline().replace('#', '').lower()).split()
cols_index = [cols_names.index(input_names[para].lower()) for para in paras]
## get values
val0 = np.loadtxt(inpath)[-1, cols_index]

# load the test results
val_list = []
for infile in infile_list:
    inpath = os.path.join(indir, infile)

    ## get column names
    with open(inpath) as tmp:
        cols_names = (tmp.readline().replace('#', '').lower()).split()
    cols_index = [cols_names.index(input_names[para].lower()) for para in paras]
    ## get values
    val_tmp = np.loadtxt(inpath)[-1, cols_index]
    val_list.append(val_tmp)

# loop over parameters
for i_para, para in enumerate(paras):

    # find the max shift
    xmin = 0
    xmax = 0
    for val_tmp in val_list:

        xval = val_tmp[i_para] - val0[i_para]
        if xval > xmax:
            xmax = xval
        if xval < xmin:
            xmin = xval

    print(f">>>>> {para}: {val0[i_para]}, shift: {xmin}, {xmax}")
