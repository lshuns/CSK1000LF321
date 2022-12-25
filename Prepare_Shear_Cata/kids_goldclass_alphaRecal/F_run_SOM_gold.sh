#!/usr/bin/bash

# -*- coding: utf-8 -*-
# @Author: ssli
# @Date:   2022-12-24 12:38:35
# @Last Modified by:   ssli
# @Last Modified time: 2022-12-25 12:18:35

### run the SOM gold selection
###### for Leiden, should run in cosmoW environment
###### for general machine, should have CosmoWrapper installed (https://github.com/lshuns/CosmoWrapper)

WORKDIR=`pwd`

bash ../utils/SOM_gold.sh ${WORKDIR}/F_config_SOM_gold.param

# ########## 1 core
# use parameter file /net/grecht/data2/ssli_files/Projects/Projects/10CSK1000/Prepare_Shear_Cata/kids_goldclass_alphaRecal/F_config_SOM_gold.param
# Constructing Spectroscopic Adapt Catalogue {
# Reading /disks/shear10/ssli/K1000CS/LF321_Inputs//KiDS_specz_PAUS_COSMOS2015.csv - Done
# MAG_GAAPadapt_Z->MAG_GAAP_Z
# MAG_GAAPadapt_Y->MAG_GAAP_Y
# MAG_GAAPadapt_J->MAG_GAAP_J
# MAG_GAAPadapt_H->MAG_GAAP_H
# MAG_GAAPadapt_Ks->MAG_GAAP_Ks
# MAGERR_GAAPadapt_Z->MAGERR_GAAP_Z
# MAGERR_GAAPadapt_Y->MAGERR_GAAP_Y
# MAGERR_GAAPadapt_J->MAGERR_GAAP_J
# MAGERR_GAAPadapt_H->MAGERR_GAAP_H
# MAGERR_GAAPadapt_Ks->MAGERR_GAAP_Ks
# MAG_GAAPadapt_u->MAG_GAAP_u
# MAG_GAAPadapt_g->MAG_GAAP_g
# MAG_GAAPadapt_r->MAG_GAAP_r
# MAG_GAAPadapt_i->MAG_GAAP_i
# Writing /disks/shear10/ssli/K1000CS/LF321_Inputs//KiDS_specz_PAUS_COSMOS2015_adapt.fits - Done
# } - Done
# Photometric Catalogue Already Exists! Skipping!
# Constructing the DIR Column Photometry Catalogue {
# > All done
# > All done
# } - Done
# Constructing the Fiducial SOM {
# Creating DIR weights with input training catalogue(s):
#    --> /disks/shear10/ssli/K1000CS/LF321_Inputs//KiDS_specz_PAUS_COSMOS2015_adapt.fits
# And corresponding input reference catalogue(s):
#    --> /disks/shear10/ssli/K1000CS/LF321_Inputs//K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A1_DIRcols.cat
# Working on Catalogues:
#      /disks/shear10/ssli/K1000CS/LF321_Inputs//KiDS_specz_PAUS_COSMOS2015_adapt.fits 
#      /disks/shear10/ssli/K1000CS/LF321_Inputs//K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A1_DIRcols.cat 
#   > Reading Training Catalogue - Done 33 sec 
#   > Reading Reference Catalogue - Done 0 sec 
#   > Generating the SOM from Training catalogue  (multiple numeric factors) {
#     -> constructing scaled data vector    -> whitening the input data - Done 1 sec  
#     -> constructing som grid (101x101) - Done 0 sec  
#     -> constructing SOM from full data vector - Done 2 hr 0 min 40 sec 
#  - Done 2 hr 0 min 43 sec  
#     -> Skipping remaining loop because '--only.som' flag was set!
#        (In the future, you could just use the kohtrain function directly)
# Warning message:
# In while (class(hdr) != "try-error") { :
#   the condition has length > 1 and only the first element will be used
# DIR weight Generation completed in 2 hr 1 min 21 sec ( 33 sec catalogue IO)

# } - Done
# Constructing the Goldclass subsets {
# Warning messages:
# 1: Modifying Blinds away from Fiducial?! 
# 2: Modifying count variable from Fiducial?! 
# 3: Modifying redshift variable from Fiducial?! 
# Warning messages:
# 1: In while (class(hdr) != "try-error") { :
#   the condition has length > 1 and only the first element will be used
# 2: In while (class(hdr) != "try-error") { :
#   the condition has length > 1 and only the first element will be used
# Warning message:
# In outputpath("_allgoldclass", phot.catalogue, output.folder) :
#   Cannot write to LDAC; preparing to write as FITS instead!
# Warning message:
# In outputpath("_allgoldclass", phot.catalogue, output.folder) :
#   Cannot write to LDAC; preparing to write as FITS instead!
# Fid: abs(mean(spec$Zbest,na.rm=T)-weighted.mean(phot$Z_B,phot$AlphaRecalC_weight,na.rm=T))<=0.61
# There were 15 warnings (use warnings() to see them)
# Warning message:
# In outputpath("_allgoldclass", phot.catalogue, output.folder) :
#   Cannot write to LDAC; preparing to write as FITS instead!

# real    73m54.557s
# user    817m22.489s
# sys 49m14.288s
# } - Done
# Constructing merged catalogue {
# Reading user blindlist: ['NONE']
# Reading user columnlist: ['Flag_SOM_Fid']
# /disks/shear10/ssli/K1000CS/LF321_Inputs//K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A1.cat
# /disks/shear10/ssli/K1000CS/LF321_Inputs//K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A1_DIRcols_allgoldclass.fits
# /disks/shear10/ssli/K1000CS/LF321_Inputs//K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A1_goldclasses.cat
# ['NONE']
# Everything is well with the world - continue
# Everything is well with the world - continue
# Looping over Columns
# adding column Flag_SOM_Fid_NONE to ldac table
# writing out file /disks/shear10/ssli/K1000CS/LF321_Inputs//K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A1_goldclasses.cat
# } - Done
# Elapsed:3:25:02.97,User=55555.173,System=3224.013,CPU=477.7%.