#!/usr/bin/bash

# -*- coding: utf-8 -*-
# @Author: ssli
# @Date:   2022-12-24 12:38:35
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-22 16:24:02

### run the SOM gold selection
###### for Leiden, should run in cosmoW environment
###### for general machine, should have CosmoWrapper installed (https://github.com/lshuns/CosmoWrapper)

WORKDIR=`pwd`

bash ../utils/SOM_gold.sh ${WORKDIR}/F_config_SOM_gold.param

# (cosmowrapper) eemmeer [125] > bash F_run_SOM_gold.sh
# use parameter file /net/grecht/data2/ssli_files/Projects/Projects/10CSK1000/A_Prepare_Shear_Cata/kids_goldclass_alphaRecal/F_config_SOM_gold.param
# Spectroscopic Adapt Catalogue Already Exists! Skipping!
# > All done
# Link the LDAC photcat {
# ‘/disks/shear10/ssli/K1000CS/LF321_Outputs//K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A1.cat’ -> ‘/disks/shear10/ssli/K1000CS/LF321_Inputs//K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A1.cat’
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
#  [1] "THELI_NAME"         "SID"                "SeqNr"             
#  [4] "MAG_AUTO"           "MAG_GAAP_u"         "MAG_GAAP_g"        
#  [7] "MAG_GAAP_r"         "MAG_GAAP_i"         "MAG_GAAP_Z"        
# [10] "MAG_GAAP_Y"         "MAG_GAAP_J"         "MAG_GAAP_H"        
# [13] "MAG_GAAP_Ks"        "Z_B"                "Z_B_MIN"           
# [16] "Z_B_MAX"            "T_B"                "AlphaRecalC_weight"
# Parse Spec Data
# Running kohparse
# Loading factor names from SOM
# Running data prediction
# Ending
# Done in 2 min 29 sec
# Parse phot Data
# Running kohparse
# Loading factor names from SOM
# Running data prediction
# Ending
# Done in 18 min 29 sec
# Warning message:
# In outputpath("_allgoldclass", phot.catalogue, output.folder) :
#   Cannot write to LDAC; preparing to write as FITS instead!
# Output File does not exist!
# Creating GoldFlag and GoldClass
# nQ≥1 only & blind NONE &  $ 880 $ &  $ 94.6 $ &  $ 98.5 $ &  $ 98.1 $ &  $ 99.4 $ & 
# nQ≥2 only & blind NONE &  $ 85.3 $ &  $ 90.2 $ &  $ 94.8 $ &  $ 94.1 $ &  $ 97.4 $ & 
# nQ≥3 only & blind NONE &  $ 83.8 $ &  $ 87.5 $ &  $ 93.1 $ &  $ 93.1 $ &  $ 97.3 $ & 
# nQ≥4 only & blind NONE &  $ 74.4 $ &  $ 81.6 $ &  $ 88.5 $ &  $ 880 $ &  $ 950 $ & 
# Warning message:
# In outputpath("_allgoldclass", phot.catalogue, output.folder) :
#   Cannot write to LDAC; preparing to write as FITS instead!


#  GoldClass & Blind & Tomographic Bin & Delta n_{\rm eff} \
# Fid: abs(mean(spec$Zbest,na.rm=T)-weighted.mean(phot$Z_B,phot$AlphaRecalC_weight,na.rm=T))<=0.61
# Fid & NONE & $ 76.8 $ &  $ 80.3 $ &  $ 860 $ &  $ 89.8 $ &  $ 950 $ & 
# There were 15 warnings (use warnings() to see them)
# Warning message:
# In outputpath("_allgoldclass", phot.catalogue, output.folder) :
#   Cannot write to LDAC; preparing to write as FITS instead!

# real    132m53.184s
# user    768m49.009s
# sys 57m55.064s
# } - Done
# Constructing merged catalogue {
# Reading user blindlist: ['NONE']
# Reading user columnlist: ['Flag_SOM_Fid']
# /disks/shear10/ssli/K1000CS/LF321_Outputs//K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A1.cat
# /disks/shear10/ssli/K1000CS/LF321_Outputs//K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A1_DIRcols_allgoldclass.fits
# /disks/shear10/ssli/K1000CS/LF321_Outputs//K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A1_goldclasses.cat
# ['NONE']
# Everything is well with the world - continue
# Everything is well with the world - continue
# Looping over Columns
# adding column Flag_SOM_Fid_NONE to ldac table
# writing out file /disks/shear10/ssli/K1000CS/LF321_Outputs//K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A1_goldclasses.cat
# } - Done
# Elapsed:2:25:44.43,User=46490.619,System=3627.901,CPU=573.1%.
