#!/usr/bin/bash

# -*- coding: utf-8 -*-
# @Author: ssli
# @Date:   2022-10-27 15:01:23
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-22 17:03:28

### run alphaRecal step2
###### require alphaRecal code, available from https://github.com/KiDS-WL/MultiBand_ImSim/tree/main/alphaRecalPlus

# usage: step2_methodD.py [-h] [--inpath INPATH] [--outDir OUTDIR]
#                         [--col_goldFlag COL_GOLDFLAG]
#                         [--col_weight COL_WEIGHT] [--col_snr COL_SNR]
#                         [--col_ZB COL_ZB] [--cols_e12 COLS_E12 COLS_E12]
#                         [--cols_psf_e12 COLS_PSF_E12 COLS_PSF_E12]
#                         [--Z_B_edges [Z_B_EDGES [Z_B_EDGES ...]]]

# step2_methodD.py: correct alpha in e1,2 with method D.

# optional arguments:
#   -h, --help            show this help message and exit
#   --inpath INPATH       the in path for the catalogue.
#                             supported formats: feather, fits and cat
#   --outDir OUTDIR       directory for the final catalogue.
#                             outpath name will be inpath_name + .A2.feather
#   --col_goldFlag COL_GOLDFLAG
#                         columns to the redshift gold class in the catalogue.
#   --col_weight COL_WEIGHT
#                         columns to the weight in the catalogue.
#   --col_snr COL_SNR     columns to the SNR in the catalogue.
#   --col_ZB COL_ZB       columns to the Z_B in the catalogue.
#   --cols_e12 COLS_E12 COLS_E12
#                         column names for e1_gal, e2_gal.
#   --cols_psf_e12 COLS_PSF_E12 COLS_PSF_E12
#                         column names for e1_psf, e2_psf.
#   --Z_B_edges [Z_B_EDGES [Z_B_EDGES ...]]
#                         edges for tomographic binning.

python /net/grecht/data2/ssli_files/Projects/Projects/8ShearBias_ImSim/MultiBand_ImSim/alphaRecalPlus/step2_methodD.py\
    --inpath /disks/shear10/ssli/K1000CS/LF321_Outputs//K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A1_goldclasses.cat\
    --outDir /disks/shear10/ssli/K1000CS/LF321_Outputs/\
    --col_goldFlag Flag_SOM_Fid_NONE\
    --col_weight AlphaRecalC_weight\
    --col_snr model_SNratio\
    --col_ZB Z_B\
    --cols_e12 autocal_e1 autocal_e2\
    --cols_psf_e12 PSF_e1 PSF_e2\
    --Z_B_edges 0.1 0.3 0.5 0.7 0.9 1.2

# (py377) eemmeer [130] > bash G_alphaRecal_step2.sh 
# number original 32067100
# number after gold selection 23401765
# number after weight selection 23401765
# number after Z_B selection 23401765
# alpha map produced in 183.9499933719635 s
# number with meaningful e after D1 23401764 fraction 0.9999999572681804
# D1 finished in 103.41310358047485 s
# number with meaningful e after D2 23401764 fraction 1.0
# D2 finished in 595.150342464447 s
# number in final cata 23401764
# final results saved to /disks/shear10/ssli/K1000CS/LF321_Outputs//K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A1_goldclasses.cat.A2.feather
# Elapsed:37:37.52,User=1232.151,System=2703.274,CPU=174.3%.