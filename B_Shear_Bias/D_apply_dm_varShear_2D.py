# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-09-29 20:23:37
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-11 17:20:38

### apply the dm from varShear to the whole sample results
####### in 2D bins (R and SNR)

import os
import re
import glob
import argparse

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

# +++++++++++++++++++++++++++++ general info

# the ZB bins
N_Zbins = 5
ZB_mins = [0.1, 0.1, 0.3, 0.5, 0.7, 0.9]
ZB_maxs = [1.2, 0.3, 0.5, 0.7, 0.9, 1.2]

# the raw results
inpath_raw_list = ['./results/m_skills_v07D7p1_LF_321_kidsPhotometry_A12_gold_Rewei_whole.csv']
inpath_raw_list += [f'./results/m_skills_v07D7p1_LF_321_kidsPhotometry_A12_gold_Rewei_bin{i_bin}.csv'
                    for i_bin in range(N_Zbins)]

# the blending-only correction
inpath_blend_list = ['./results/dm_forBlending_dz0p1_whole.csv']
inpath_blend_list += [f'./results/dm_forBlending_dz0p1_bin{i_bin}.csv' for i_bin in range(N_Zbins)]

# the blending fraction
inpath_fraction_list = ['./results/blending_fraction_whole.csv']
inpath_fraction_list += [f'./results/blending_fraction_bin{i_bin}.csv' for i_bin in range(N_Zbins)]

# where to save
outpath = './results/m_skills_v07D7p1_LF_321_kidsPhotometry_A12_gold_Rewei_varCorr_dz0p1_2D.csv'

# +++++++++++++++++++++++++++ correct the m 

res_df = pd.DataFrame(-999, 
                index = np.arange(N_Zbins+1), 
                columns = ['m1', 'm2', 'c1', 'c2',
                        'm1_err', 'm2_err', 'c1_err', 'c2_err',
                        'Z_B_min', 'Z_B_max', 'mean_blendFraction'])

for i_path, inpath_raw in enumerate(inpath_raw_list):

    # the range info
    res_df.loc[i_path, 'Z_B_min'] = ZB_mins[i_path]
    res_df.loc[i_path, 'Z_B_max'] = ZB_maxs[i_path]

    inpath_blend = inpath_blend_list[i_path]
    inpath_fraction = inpath_fraction_list[i_path]

    cata_raw = pd.read_csv(inpath_raw)
    cata_blend = pd.read_csv(inpath_blend)
    cata_fraction = pd.read_csv(inpath_fraction)

    # sort for easy combine
    cata_raw.sort_values(by=['bin1_id', 'bin2_id'], inplace=True, ignore_index=True)
    cata_blend.sort_values(by=['bin1_id', 'bin2_id'], inplace=True, ignore_index=True)
    cata_fraction.sort_values(by=['bin1_id', 'bin2_id'], inplace=True, ignore_index=True)

    # results after correction
    ratio_blend = cata_fraction['WeiBlend'].values/cata_fraction['WeiTot'].values
    for col_m in ['m1', 'm2']:
        res_df.loc[i_path, col_m] = np.average(cata_raw[col_m].values \
                                        + cata_blend[col_m].values * ratio_blend,
                                        weights=cata_raw['dataWei'].values)
        res_df.loc[i_path, f'{col_m}_err'] = np.sqrt(
                                                np.sum(
                                                    np.square(cata_raw['dataWei'].values \
                                                        * np.hypot(cata_raw[f'{col_m}_err'].values, 
                                                            cata_blend[f'{col_m}_err'].values * ratio_blend)
                                                        )
                                                    )
                                                ) / np.sum(cata_raw['dataWei'].values)
    ## c is not changed
    for col_c in ['c1', 'c2']:
        res_df.loc[i_path, col_c] = np.average(cata_raw[col_c].values,
                                        weights=cata_raw['dataWei'].values)
        res_df.loc[i_path, f'{col_c}_err'] = np.sqrt(
                                                np.sum(
                                                    np.square(cata_raw['dataWei'].values \
                                                        * cata_raw[f'{col_c}_err'].values
                                                        )
                                                    )
                                                ) / np.sum(cata_raw['dataWei'].values)
    ## the average blending fraction
    res_df.loc[i_path, 'mean_blendFraction'] = np.average(ratio_blend,
                                        weights=cata_raw['simWei'].values)

# save
print(res_df)
res_df.to_csv(outpath, index=False)
print('saved to', outpath)