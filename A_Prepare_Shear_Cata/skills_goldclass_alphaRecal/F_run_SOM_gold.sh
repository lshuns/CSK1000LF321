#!/usr/bin/bash

# -*- coding: utf-8 -*-
# @Author: ssli
# @Date:   2022-12-24 12:38:35
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-22 16:25:51

### run the SOM gold selection
###### for Leiden, should run in cosmoW environment
###### for general machine, should have CosmoWrapper installed (https://github.com/lshuns/CosmoWrapper)

WORKDIR=`pwd`

bash ../utils/SOM_gold.sh ${WORKDIR}/F_config_SOM_gold.param

# (cosmowrapper) markermeer [110] > bash F_run_SOM_gold.sh 
# use parameter file /net/grecht/data2/ssli_files/Projects/Projects/10CSK1000/A_Prepare_Shear_Cata/skills_goldclass_alphaRecal/F_config_SOM_gold.param
# Spectroscopic Adapt Catalogue Already Exists! Skipping!
# > All done
# Link the LDAC photcat {
# ‘/disks/shear10/ssli/K1000CS//skills_v07D7p1_Outputs//skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_newCut_A1_WeiCut.cat’ -> ‘/disks/shear10/ssli/K1000CS//skills_v07D7p1_Inputs//skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_newCut_A1_WeiCut.cat’
# } - Done
# Constructing the DIR Column Photometry Catalogue {
# > All done
# > All done
# } - Done
# Fiducial SOM Already Exists! Skipping!
# Constructing the Goldclass subsets {
# WARNING: Modifying Blinds away from Fiducial?!
# WARNING: Modifying count variable from Fiducial?!
# WARNING: Modifying redshift from Fiducial?!
# Warning messages:
# 1: Modifying Blinds away from Fiducial?! 
# 2: Modifying count variable from Fiducial?! 
# 3: Modifying redshift variable from Fiducial?! 
# Load SOM Data
# Read Phot Data
# Warning messages:
# 1: In while (class(hdr) != "try-error") { :
#   the condition has length > 1 and only the first element will be used
# 2: In while (class(hdr) != "try-error") { :
#   the condition has length > 1 and only the first element will be used
#  [1] "MAG_AUTO"           "MAG_GAAP_u"         "MAG_GAAP_g"        
#  [4] "MAG_GAAP_r"         "MAG_GAAP_i"         "MAG_GAAP_Z"        
#  [7] "MAG_GAAP_Y"         "MAG_GAAP_J"         "MAG_GAAP_H"        
# [10] "MAG_GAAP_Ks"        "Z_B"                "T_B"               
# [13] "AlphaRecalC_weight"
# Parse Spec Data
# Running kohparse
# Loading factor names from SOM
# Running data prediction
# Ending
# Done in 2 min 26 sec
# Parse phot Data
# Running kohparse
# Loading factor names from SOM
# Running data prediction
# Ending
# Done in 22 min 21 sec
# Warning message:
# In outputpath("_allgoldclass", phot.catalogue, output.folder) :
#   Cannot write to LDAC; preparing to write as FITS instead!
# Output File does not exist!
# Creating GoldFlag and GoldClass
# nQ≥1 only & blind NONE &  $ 94.7 $ &  $ 96.9 $ &  $ 98.6 $ &  $ 98.5 $ &  $ 99.3 $ & 
# nQ≥2 only & blind NONE &  $ 91.8 $ &  $ 94.2 $ &  $ 95.1 $ &  $ 960 $ &  $ 97.5 $ & 
# nQ≥3 only & blind NONE &  $ 90.5 $ &  $ 920 $ &  $ 93.5 $ &  $ 94.8 $ &  $ 97.3 $ & 
# nQ≥4 only & blind NONE &  $ 83.2 $ &  $ 86.2 $ &  $ 88.4 $ &  $ 90.5 $ &  $ 94.6 $ & 
# Warning message:
# In outputpath("_allgoldclass", phot.catalogue, output.folder) :
#   Cannot write to LDAC; preparing to write as FITS instead!


#  GoldClass & Blind & Tomographic Bin & Delta n_{\rm eff} \
# Fid: abs(mean(spec$Zbest,na.rm=T)-weighted.mean(phot$Z_B,phot$AlphaRecalC_weight,na.rm=T))<=0.61
# Fid & NONE & $ 79.4 $ &  $ 85.8 $ &  $ 85.2 $ &  $ 92.6 $ &  $ 94.7 $ & 
# There were 15 warnings (use warnings() to see them)
# Warning message:
# In outputpath("_allgoldclass", phot.catalogue, output.folder) :
#   Cannot write to LDAC; preparing to write as FITS instead!

# real    122m43.617s
# user    876m47.924s
# sys 64m40.365s
# } - Done
# Constructing merged catalogue {
# Reading user blindlist: ['NONE']
# Reading user columnlist: ['Flag_SOM_Fid']
# /disks/shear10/ssli/K1000CS//skills_v07D7p1_Outputs//skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_newCut_A1_WeiCut.cat
# /disks/shear10/ssli/K1000CS//skills_v07D7p1_Outputs//skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_newCut_A1_WeiCut_DIRcols_allgoldclass.fits
# /disks/shear10/ssli/K1000CS//skills_v07D7p1_Outputs//skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_newCut_A1_WeiCut_goldclasses.cat
# ['NONE']
# Everything is well with the world - continue
# Everything is well with the world - continue
# Looping over Columns
# adding column Flag_SOM_Fid_NONE to ldac table
# writing out file /disks/shear10/ssli/K1000CS//skills_v07D7p1_Outputs//skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_newCut_A1_WeiCut_goldclasses.cat
# } - Done
# Elapsed:2:13:09.37,User=52947.507,System=3987.045,CPU=712.6%.
