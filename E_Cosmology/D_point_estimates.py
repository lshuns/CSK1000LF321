# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-05-22 11:53:59
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-22 12:30:35

### get the individual parameter constraints from the MCMC chain
###### three statistics: PJ-HPD, Mean, Max
###### PJ-HPD reference: Section 6.4 of Joachimi et al. (2021)
###### Mean from postprocess within CosmoSIS  
###### Max from ChainConsumer with statistics='max' and kde=1.0

import os
import sys
import pandas as pd 
import numpy as np 

from chainconsumer import ChainConsumer

# code from kcap for calculating PJ-HPD
sys.path.append('./external')
from stat_tools import find_projected_joint_HPDI

# >>>>>>>>>>>>>>>>>>>> I/O

run_names = ["multinest_K1000_like", "polychord_nIA", "multinest_nIA", "multinest_nIA_noCut", "multinest_wIA", "multinest_wIA_eta", "multinest_nIA_noCut", "multinest_wIA_noCut"]
LABELs = ['MA', 'JLvdB'] + run_names
## the MCMC chain files
inpath_chain_list = ['/disks/shear10/ssli/K1000CS/K1000_info/KiDS1000_cosmis_shear_data_release/chains_and_config_files/main_chains_iterative_covariance/cosebis/chain/output_multinest_C.txt',
'/disks/shear10/ssli/K1000CS/K1000_info/KiDS1000_vdB22_cosmic_shear_data_release/multinest/Fid_output_multinest_C.txt']
for run_name in run_names:
    if 'polychord' in run_name:
        prefix = 'chain'
    else:
        prefix = 'output'
    inpath_chain_list.append(f'/disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/{run_name}/{prefix}_{run_name}.txt')
## the local minimisation file
inpath_max_list = ['/disks/shear10/ssli/K1000CS/K1000_info/KiDS1000_cosmis_shear_data_release/chains_and_config_files/main_chains_iterative_covariance/cosebis/chain/maxpost_multinest_start_C.txt',
'/disks/shear10/ssli/K1000CS/K1000_info/KiDS1000_vdB22_cosmic_shear_data_release/maxlike/Fid_output_maxlike_C.txt']\
+ [f'/disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/{run_name}/maxpost_{run_name}_start.txt' for run_name in run_names]
## the CosmoSIS summary
indir_summary_list = ['/disks/shear10/ssli/K1000CS/K1000_info/KiDS1000_cosmis_shear_data_release/chains_and_config_files/main_chains_iterative_covariance/cosebis/chain/output_multinest_C_summary',
'/disks/shear10/ssli/K1000CS/K1000_info/KiDS1000_vdB22_cosmic_shear_data_release/multinest/Fid_output_multinest_C_summary']\
+ [f'/disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/{run_name}/cosmosis_summary' for run_name in run_names]

## which parameters are required
paras = ['s8', 'omega_m', 'sigma_8', 'a_ia']

## the level of credible region
alpha = 0.683

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

# initialise chainconsumer
c = ChainConsumer()

# get parameter names in chain
para_names = [input_names[para].lower() for para in paras]

# collect values
MAP_list = []
CI_list = []
mean_list = []
lerr68_list = []
uerr68_list = []
for i_inpath, label in enumerate(LABELs):
    # print(">>>>> for ", label)
    inpath_chain = inpath_chain_list[i_inpath]
    inpath_max = inpath_max_list[i_inpath]
    indir_summary = indir_summary_list[i_inpath]

    # load the summary from cosmosis
    mean_l_u = []
    for sum_file in ['means.txt', 'lerr68.txt', 'uerr68.txt']:
        tmp_list = []
        with open(os.path.join(indir_summary, sum_file), 'r') as tmp:
            lines = tmp.readlines() 
        for para_name in para_names:
            for index, line in enumerate(lines):
                para_name0 = line.split()[0].lower()
                if para_name == para_name0:
                    tmp_list.append(float(line.split()[1]))
                    break
        del lines
        mean_l_u.append(tmp_list)
        del tmp_list
    mean_list.append(mean_l_u[0])
    lerr68_list.append(mean_l_u[1])
    uerr68_list.append(mean_l_u[2])
    del mean_l_u

    # load the chain
    ## the col names
    with open(inpath_chain) as tmp:
        cols_chain = (tmp.readline().replace('#', '').lower()).split()
    ## the values
    vals_chain0 = np.loadtxt(inpath_chain)
    ## get the column index for the required parameters
    cols_index_chain = [cols_chain.index(input_names[para].lower()) for para in paras]
    ## weights and posterior
    log_posterior = vals_chain0[:, cols_chain.index('post')]
    weights = vals_chain0[:, cols_chain.index('weight')]
    ## get the values for the required parameters
    vals_chain = [vals_chain0[:, col_index] for col_index in cols_index_chain]

    # add to chainconsumer
    c.add_chain(vals_chain0[:, cols_index_chain], weights=weights, posterior=log_posterior,
            parameters=paras, name=label)
    del vals_chain0, cols_index_chain

    # load MAP if given
    #### assuming the last row has the max post
    if inpath_max is not None:
        with open(inpath_max) as tmp:
            cols_max = (tmp.readline().replace('#', '').lower()).split()
        vals_max = np.loadtxt(inpath_max)[-1]
        cols_index_max = [cols_max.index(input_names[para].lower()) for para in paras]
        # get the MAP estimate
        vals_max = [vals_max[col_index] for col_index in cols_index_max]
        del cols_index_max
    else:
        vals_max = None

    CI_list_tmp = []
    MAP_list_tmp = []
    for i_para, para in enumerate(paras):

        samples = vals_chain[i_para]
        if vals_max is not None:
            MAP = vals_max[i_para]
        else:
            MAP = None

        # get the PJ-HPD
        CI, MAP_chain, coverage_1d, n_sample = find_projected_joint_HPDI(samples, weights=weights, coverage_1d_threshold=alpha, 
                                      MAP=MAP,
                                      sort_idx=None, log_posterior=log_posterior, 
                                      method="interpolate", twosided=True, verbose=False, strict=False,
                                      return_map=True, 
                                      return_coverage_1d=True, return_coverage_nd=False, 
                                      return_n_sample=True)

        CI_list_tmp.append(CI)
        if MAP is not None:
            ### the one reported in Asgari 2021
            MAP_list_tmp.append(MAP)
        else:
            MAP_list_tmp.append(MAP_chain)
    CI_list.append(CI_list_tmp)
    MAP_list.append(MAP_list_tmp)

# the max
c.configure(kde=1., statistics="max")
summary_max_list = c.analysis.get_summary()

# print results
for i_inpath, label in enumerate(LABELs):
    print(">>>>>>>>>>>>>>>>> for", label)

    for i_para, para in enumerate(paras):
        print("+++ parameter:", para)

        # the PJ-HPD
        res_tmp = [ MAP_list[i_inpath][i_para], 
                    CI_list[i_inpath][i_para][0] - MAP_list[i_inpath][i_para],
                    CI_list[i_inpath][i_para][1] - MAP_list[i_inpath][i_para], 
                    (CI_list[i_inpath][i_para][1] - CI_list[i_inpath][i_para][0])/2.]
        print('====== MAP', ', '.join('{:.3f}'.format(res) for res in res_tmp))
        print('========== latex', f'${res_tmp[0]:.3f}_{{{res_tmp[1]:.3f}}}^{{+{res_tmp[2]:.3f}}}$')

        # the mean
        res_tmp = [mean_list[i_inpath][i_para],
                    lerr68_list[i_inpath][i_para] - mean_list[i_inpath][i_para],
                    uerr68_list[i_inpath][i_para] - mean_list[i_inpath][i_para], 
                    (uerr68_list[i_inpath][i_para] - lerr68_list[i_inpath][i_para])/2.]
        print('====== Mean', ', '.join('{:.3f}'.format(res) for res in res_tmp))
        print('========== latex', f'${res_tmp[0]:.3f}_{{{res_tmp[1]:.3f}}}^{{+{res_tmp[2]:.3f}}}$')

        # the max
        res_tmp = [summary_max_list[i_inpath][para][1],
                    summary_max_list[i_inpath][para][0] - summary_max_list[i_inpath][para][1],
                    summary_max_list[i_inpath][para][2] - summary_max_list[i_inpath][para][1], 
                    (summary_max_list[i_inpath][para][2] - summary_max_list[i_inpath][para][0])/2.]
        print('====== Max', ', '.join('{:.3f}'.format(res) for res in res_tmp))
        print('========== latex', f'${res_tmp[0]:.3f}_{{{res_tmp[1]:.3f}}}^{{+{res_tmp[2]:.3f}}}$')
