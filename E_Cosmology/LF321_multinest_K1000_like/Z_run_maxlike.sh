#!/usr/bin/bash

# @Author: lshuns
# @Date:   2023-05-01 10:39:51
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-05 18:49:56

### run local minimisation with Nelder-Mead method

# activate kcap
source /disks/shear10/ssli/K1000CS/my_kcap/env/bin/activate
source cosmosis-configure

# run Marika's code
python ../maxlike/maxlike_cosmosis.py -i ./pipeline.ini \
  -m /disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/multinest/output_multinest.txt\
  -o /disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/multinest/maxpost_multinest_start.txt\
  --max_post -s \
  --maxiter 3000

# Nelder-Mead results:
# cosmological_parameters--omch2
# cosmological_parameters--omch2 0.09279699171927089 0.08938588099736283
# cosmological_parameters--ombh2
# cosmological_parameters--ombh2 0.022187833272010758 0.023216973052777104
# cosmological_parameters--h0
# cosmological_parameters--h0 0.6443098761165579 0.6441012109516991
# cosmological_parameters--n_s
# cosmological_parameters--n_s 1.082877664835508 1.0779405256873555
# cosmological_parameters--s_8_input
# cosmological_parameters--s_8_input 0.7689312800826512 0.7726293931573006
# halo_model_parameters--a
# halo_model_parameters--a 2.2144091609867997 2.432848222058806
# intrinsic_alignment_parameters--a
# intrinsic_alignment_parameters--a 0.25162078008178757 0.22731188172172345
# nofz_shifts--uncorr_bias_1
# nofz_shifts--uncorr_bias_1 0.08989152341300599 -0.03805935812465577
# nofz_shifts--uncorr_bias_2
# nofz_shifts--uncorr_bias_2 0.5789918635584961 0.40937080109301865
# nofz_shifts--uncorr_bias_3
# nofz_shifts--uncorr_bias_3 -1.694256301157302 -1.7642981230517305
# nofz_shifts--uncorr_bias_4
# nofz_shifts--uncorr_bias_4 -1.6026987501824141 -1.958267008516259
# nofz_shifts--uncorr_bias_5
# nofz_shifts--uncorr_bias_5 1.466145998787967 1.0670643120113257
# min chi2= 61.409841886811854
# nofz_shifts
# Elapsed:2:35:24.19,User=101855.514,System=97.565,CPU=1093.4%.
