# -*- coding: utf-8 -*-
# @Author: ssli
# @Date:   2022-10-27 15:56:26
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-22 22:11:21

### measure correlation function from SKiLLS

import os
import shutil
import pathlib

import numpy as np
import pandas as pd 

from CorrFunc import CalCorrFunc

# >>>>>> I/O

inpath_list = ['/disks/shear10/ssli/K1000CS/skills_v07D7p1_Outputs/skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_newCut_A1_WeiCut_goldclasses.cat.A2.feather',
            '/disks/shear10/ssli/K1000CS/skills_v07D7p1_Outputs/skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_newCut_A1_WeiCut_goldclasses.cat.A2.feather.shear_m283m283.feather',
            '/disks/shear10/ssli/K1000CS/skills_v07D7p1_Outputs/skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_newCut_A1_WeiCut_goldclasses.cat.A2.feather.shear_m283p283.feather',
            '/disks/shear10/ssli/K1000CS/skills_v07D7p1_Outputs/skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_newCut_A1_WeiCut_goldclasses.cat.A2.feather.shear_p283m283.feather',
            '/disks/shear10/ssli/K1000CS/skills_v07D7p1_Outputs/skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_newCut_A1_WeiCut_goldclasses.cat.A2.feather.shear_p283p283.feather']
outfolder_list = ['CF_A12_goldclasses_whole', 
                    'CF_A12_goldclasses_m283m283',
                    'CF_A12_goldclasses_m283p283',
                    'CF_A12_goldclasses_p283m283',
                    'CF_A12_goldclasses_p283p283'
                    ]

bias_file = '../../B_Shear_Bias/results/m_skills_v07D7p1_LF_321_kidsPhotometry_A12_gold.csv'

# weight column
col_gold = 'Flag_SOM_Fid_NONE'
col_wei = 'AlphaRecalC_weight'
col_e1 = 'AlphaRecalD2_e1'
col_e2 = 'AlphaRecalD2_e2'

# >>>>>> general setups

# reshift bins
ZB_bins = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2]

# parameters for treecorr
theta_nbins = 9
theta_min = 0.5
theta_max = 300.
theta_unit = "arcmin"
theta_bin_slop = 0.05
nthr = 64

# used columns
cols_used = [col_gold, col_wei, col_e1, col_e2, 
                'Z_B', 'X_WORLD', 'Y_WORLD', 
                ]
cols_rename = {col_wei: 'shape_weight', 
                col_e1: 'e1',
                col_e2: 'e2',
                'X_WORLD': 'RA',
                'Y_WORLD': 'DEC'}

# >>>>>> get the c term
tmp = pd.read_csv(bias_file)
c1_list = tmp.loc[1:, 'c1'].values
print('c1', c1_list)
c2_list = tmp.loc[1:, 'c2'].values
print('c2', c2_list)

# >>>>>> workhorse

for i_path, inpath in enumerate(inpath_list):

    print('>>> working on', inpath)

    outfolder = outfolder_list[i_path]
    outdir = os.path.join(os.path.dirname(inpath), outfolder)
    print('results will saved to', outdir)

    # 0. clean previous results everytime
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.mkdir(outdir)

    # 1. load the catalogue
    cata = pd.read_feather(inpath)
    print('number original', len(cata))
    ## used columns
    cata = cata[cols_used]
    ## rename
    cata.rename(columns=cols_rename, inplace=True)
    ## only within the tomo bins
    cata = cata[(cata['shape_weight']>0)]
    cata.reset_index(drop=True, inplace=True)
    print('number wei>0', len(cata))
    ## only gold class
    cata = cata[(cata[col_gold]>0)]
    cata.reset_index(drop=True, inplace=True)
    print('number goldclass', len(cata))
    ## only within the tomo bins
    cata = cata[(cata['Z_B']>ZB_bins[0])&(cata['Z_B']<=ZB_bins[-1])]
    cata.reset_index(drop=True, inplace=True)
    print('number within range', len(cata))

    # 2. bin to redshift bins
    cata.loc[:, 'i_ZBbin'] = pd.cut(cata['Z_B'].values, bins=ZB_bins, 
                                    right=True, labels=False)
    cata = cata.groupby(by='i_ZBbin')

    # 3. calculate correlation function
    Nzbins = cata.ngroups
    # print('>>> total number of redshift bins', Nzbins)
    cata_final = []
    for i_ZBbin in range(Nzbins):

        cata1 = cata.get_group(i_ZBbin)
        c12_cata1 = [c1_list[i_ZBbin], c2_list[i_ZBbin]]
        # print('+++ ZB', i_ZBbin, np.min(cata1['Z_B']), np.max(cata1['Z_B']))

        for j_ZBbin in range(i_ZBbin, Nzbins):

            # where to save
            outpath = os.path.join(outdir, f'xi_tomo{i_ZBbin}{j_ZBbin}.txt')

            # cross correlation
            if i_ZBbin != j_ZBbin:
                cata2 = cata.get_group(j_ZBbin)
                c12_cata2 = [c1_list[j_ZBbin], c2_list[j_ZBbin]]

            # auto correlation
            else:
                cata2 = None
                c12_cata2 = None

            # run
            if 'var_method' in locals():
                outpath_cov = os.path.join(outdir, f'cov_tomo{i_ZBbin}{j_ZBbin}.txt')
                CalCorrFunc(outpath,
                        cata1, cata2=cata2, 
                        theta_Nbin=theta_nbins, theta_bin_slop=theta_bin_slop, 
                        theta_min=theta_min, theta_max=theta_max, theta_unit=theta_unit, 
                        theta_bin_type='Log',
                        num_threads=nthr, 
                        c12_cata1=c12_cata1, 
                        c12_cata2=c12_cata2,
                        var_method=var_method,
                        Npatch=Npatch,
                        outpath_cov=outpath_cov)
            else:
                CalCorrFunc(outpath,
                        cata1, cata2=cata2, 
                        theta_Nbin=theta_nbins, theta_bin_slop=theta_bin_slop, 
                        theta_min=theta_min, theta_max=theta_max, theta_unit=theta_unit, 
                        theta_bin_type='Log',
                        num_threads=nthr, 
                        c12_cata1=c12_cata1, 
                        c12_cata2=c12_cata2)

            # extract used values
            cata_tmp = np.loadtxt(outpath)
            Nrows = len(cata_tmp[:, 0])
            cata_tmp = pd.DataFrame({'ito': (i_ZBbin * np.ones(Nrows)).astype(int),
                                    'jto': (j_ZBbin * np.ones(Nrows)).astype(int),
                                    'theta': cata_tmp[:, 1],
                                    'xi_p': cata_tmp[:, 3],
                                    'xi_m': cata_tmp[:, 4],
                                    'xierr_p': cata_tmp[:, 7],
                                    'xierr_m': cata_tmp[:, 8]
                                    })
            cata_final.append(cata_tmp)
    # combine and save
    cata_final = pd.concat(cata_final)
    outpath = os.path.join(outdir, f'xi_combined.csv')
    cata_final.to_csv(outpath, index=False)
    print('combined cata saved to', outpath)