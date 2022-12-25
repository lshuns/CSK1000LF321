#!/usr/bin/bash

# -*- coding: utf-8 -*-
# @Author: ssli
# @Date:   2022-10-17 10:54:26
# @Last Modified by:   ssli
# @Last Modified time: 2022-12-24 13:06:50

### run alphaRecal step1
###### require alphaRecal code, available from https://github.com/KiDS-WL/MultiBand_ImSim/tree/main/alphaRecalPlus

# ### running script
# usage: step1_methodC.py [-h] [--inpath INPATH] [--outDir OUTDIR]
#                         [--col_weight COL_WEIGHT] [--col_var COL_VAR]
#                         [--col_snr COL_SNR] [--cols_e12 COLS_E12 COLS_E12]
#                         [--cols_e12_corr COLS_E12_CORR COLS_E12_CORR]
#                         [--cols_e12_raw COLS_E12_RAW COLS_E12_RAW]
#                         [--cols_psf_e12 COLS_PSF_E12 COLS_PSF_E12]

# >>>>>>>>>>>>>>>>> KiDS 321
python /net/grecht/data2/ssli_files/Projects/Projects/8ShearBias_ImSim/MultiBand_ImSim/alphaRecalPlus/step1_methodC.py\
    --inpath /disks/shear16/ssli/KiDS/K1000_LF_321_mosaic/KIDS_1006tiles_r_SDSS.V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_selec_noWeiCut.feather\
    --outDir /disks/shear16/ssli/KiDS/K1000_LF_321_mosaic\
    --col_weight weight\
    --col_var 2D_measurement_variance\
    --col_snr model_SNratio\
    --cols_e12 autocal_e1 autocal_e2\
    --cols_e12_raw raw_e1 raw_e2\
    --cols_psf_e12 PSF_e1 PSF_e2

# number total 38,026,162
# final results saved to /disks/shear16/ssli/KiDS/K1000_LF_321_mosaic/KIDS_1006tiles_r_SDSS.V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_selec_noWeiCut_A1.feather
# Elapsed:37:08.25,User=1350.235,System=3378.901,CPU=212.2%.