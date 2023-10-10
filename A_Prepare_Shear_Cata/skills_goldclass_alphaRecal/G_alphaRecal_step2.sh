#!/usr/bin/bash

# @Author: ssli
# @Date:   2022-10-25 15:17:05
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-22 16:49:22

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
    --inpath /disks/shear10/ssli/K1000CS//skills_v07D7p1_Outputs//skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_newCut_A1_WeiCut_goldclasses.cat\
    --outDir /disks/shear10/ssli/K1000CS/skills_v07D7p1_Outputs/\
    --col_goldFlag Flag_SOM_Fid_NONE\
    --col_weight AlphaRecalC_weight\
    --col_snr SNR_LF_r\
    --col_ZB Z_B\
    --cols_e12 e1_LF_r e2_LF_r\
    --cols_psf_e12 psf_e1_LF_r psf_e2_LF_r\
    --Z_B_edges 0.1 0.3 0.5 0.7 0.9 1.2

# (py377) markermeer [112] > bash G_alphaRecal_step2.sh 
# number original 39103150
# number after gold selection 28885362
# number after weight selection 28885362
# number after Z_B selection 28885362
# alpha map produced in 149.55487704277039 s
# number with meaningful e after D1 28885361 fraction 0.9999999653803888
# D1 finished in 83.27967381477356 s
# number with meaningful e after D2 28885361 fraction 1.0
# D2 finished in 334.808908700943 s
# number in final cata 28885361
# final results saved to /disks/shear10/ssli/K1000CS//skills_v07D7p1_Outputs//skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_newCut_A1_WeiCut_goldclasses.cat.A2.feather
# Elapsed:21:33.59,User=1006.727,System=1933.893,CPU=227.3%.
