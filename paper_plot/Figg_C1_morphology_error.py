# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-04-02 17:20:20
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-09 17:27:15

### the relative uncertainties in galaxy morphology


import os

import pandas as pd 
import numpy as np 

import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True

import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator, LogLocator, NullFormatter, NullLocator

########## general info

# plot properties
# outpath = 'show'
outpath = './plots/morphology_error.pdf'
FIGSIZE = [12, 12]
font_size = 16
usetex = True
MS = 4
# font size
plt.rc('font', size=font_size)
# tex
plt.rcParams["text.usetex"] = usetex
plt.rcParams['font.family'] = 'serif'

if outpath != 'show':
    backend_orig = plt.get_backend()
    plt.switch_backend("agg")
# how many sub-plots
fig, axs = plt.subplots(3, 3, sharex=False, sharey=False, figsize=FIGSIZE)

# morphology columns
skills_cols = ['Re_arcsec', 'BA', 'shape/sersic_n']
cosmos_cols = ['RE_GALFIT_HI', 'BA_GALFIT_HI', 'N_GALFIT_HI']
err_cols = ['REERR_GALFIT_HI', 'BAERR_GALFIT_HI', 'NERR_GALFIT_HI']
dm_label_names = ['size', 'q', 'n']

# binning info for morphology columns
nbins = 40
loc_legend = None
XLABEL_h_list = ['half-light radius (arcsec)', 'axis ratio', r'S\'ersic index']
XRANGE_h_list = [[0.05, 2], [0.05, 0.9], [0.6, 5.9]]
xlog_list = [True, False, False]

# magnitude columns
cosmos_mag = 'r_mag_auto'
skills_mag = 'r_SDSS_apparent_corr'

########## COSMOS data catalogue
inpath = '/disks/shear10/ssli/ImSim/input/COSMOS_cata/cosmos_shape_z_uBVriZYJHKs.feather'
cosmos_cata = pd.read_feather(inpath)
print('COSMOS ori', len(cosmos_cata))
### select
###### 0. discard too small or too big galaxies
mask_re = (cosmos_cata['RE_GALFIT_HI']>=1e-2) & (cosmos_cata['RE_GALFIT_HI']<=10.)
###### 1. good shape
mask_galfit = (cosmos_cata['FLAG_GALFIT_HI']==0)
###### 2. has magnitude
mask_mag = (cosmos_cata['r_mag_auto']>0)
### apply
cosmos_cata = cosmos_cata[mask_galfit & mask_re & mask_mag]
del mask_galfit, mask_re
cosmos_cata.reset_index(drop=True, inplace=True)
print('COSMOS selected', len(cosmos_cata))
### selected used columns
cosmos_cata = cosmos_cata[cosmos_cols+err_cols+[cosmos_mag]]

########## SKiLLS mock catalogue
inpath = '/disks/shear10/ssli/ImSim/input/SURFS_cata/skills_v07Ds_input_part0_shifted.feather'
skills_cata = pd.read_feather(inpath)

########## dm files
indir_tmp = '/net/grecht/data2/ssli_files/Projects/Projects/8ShearBias_ImSim/MultiBand_ImSim/sensitivity_test/galaxy/results'
inpath_dm_list = [[os.path.join(indir_tmp, f'dm_ZBbins_skills_v07D7_{label}U_nogold_reweighted.csv'),  
                    os.path.join(indir_tmp, f'dm_ZBbins_skills_v07D7_{label}D_nogold_reweighted.csv')]
                    for label in dm_label_names]
COLORs_dm = ['darkred', 'darkblue']
LABELs_dm = ['shift up', 'shift down']

########## 1. the relation between the mag and relative error
i_row = 0
XLABEL = r'$r$-band magnitude'
YLABEL = 'Relative uncertainties'

XRANGE = [22.2, 25.4]
YRANGE = [0.01, 0.13]

for i_col, cosmos_col in enumerate(cosmos_cols):

    err_col = err_cols[i_col]

    ax = axs[i_row, i_col]

    ## data
    test_df = pd.DataFrame({'xval': cosmos_cata[cosmos_mag].values,
                            'yval': cosmos_cata[err_col].values/cosmos_cata[cosmos_col].values})
    # group 
    Ngroup = 10
    test_df.loc[:, 'bin'] = pd.qcut(test_df['xval'].values, Ngroup, labels=False)
    # get the median
    test_df_median = test_df.groupby('bin').median()
    ax.errorbar(test_df_median['xval'].values, test_df_median['yval'].values, 
                color='k', linestyle='--', marker='o', markersize=MS)

    # the labels
    ax.set_xlabel(XLABEL)
    if i_col == 0:
        ax.set_ylabel(YLABEL)

    # some general setting
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(XRANGE[0], XRANGE[1])
    ax.set_ylim(YRANGE[0], YRANGE[1])

########## 2. the distribution 
i_row = 1

YLABEL = 'Probability density'

# 0, 1, 2
COLORs = ['k', 'darkred', 'darkblue']
LINEs = ['solid', 'solid', 'solid']
LWs = [1, 1, 1]
LABELs = ['fiducial', 'shift up', 'shift down']

for i_col, skills_col in enumerate(skills_cols):

    ax = axs[i_row, i_col]

    xlog = xlog_list[i_col]
    XRANGE = XRANGE_h_list[i_col]
    XLABEL = XLABEL_h_list[i_col]

    # the fiducial
    para = skills_cata[skills_col].values
    if xlog:
        logbins = np.logspace(np.log10(XRANGE[0]), np.log10(XRANGE[1]), nbins)
        ax.hist(x=para, 
            bins=logbins, 
            range=XRANGE,
            density=True, 
            color=COLORs[0], 
            label=LABELs[0], 
            alpha=0.3, 
            linestyle=LINEs[0],
            linewidth=LWs[0])
    else:
        ax.hist(x=para, 
            bins=nbins, 
            range=XRANGE,
            density=True, 
            color=COLORs[0], 
            label=LABELs[0], 
            alpha=0.3, 
            linestyle=LINEs[0],
            linewidth=LWs[0])

    # the up
    para = skills_cata[f'{skills_col}_U'].values
    if xlog:
        logbins = np.logspace(np.log10(XRANGE[0]), np.log10(XRANGE[1]), nbins)
        ax.hist(x=para, 
            bins=logbins, 
            range=XRANGE,
            density=True, 
            color=COLORs[1], 
            label=LABELs[1], 
            histtype='step', 
            linestyle=LINEs[1],
            linewidth=LWs[1])
    else:
        ax.hist(x=para, 
            bins=nbins, 
            range=XRANGE,
            density=True, 
            color=COLORs[1], 
            label=LABELs[1], 
            histtype='step', 
            linestyle=LINEs[1],
            linewidth=LWs[1])

    # the down
    para = skills_cata[f'{skills_col}_D'].values
    if xlog:
        logbins = np.logspace(np.log10(XRANGE[0]), np.log10(XRANGE[1]), nbins)
        ax.hist(x=para, 
            bins=logbins, 
            range=XRANGE,
            density=True, 
            color=COLORs[2], 
            label=LABELs[2], 
            histtype='step', 
            linestyle=LINEs[2],
            linewidth=LWs[2])
    else:
        ax.hist(x=para, 
            bins=nbins, 
            range=XRANGE,
            density=True, 
            color=COLORs[2], 
            label=LABELs[2], 
            histtype='step', 
            linestyle=LINEs[2],
            linewidth=LWs[2])

    if xlog:
        ax.set_xscale('log')

    # the labels
    ax.set_xlabel(XLABEL)
    if i_col == 0:
        ax.set_ylabel(YLABEL)

    ax.set_xlim(XRANGE[0], XRANGE[1])

    if i_col==2:
        ax.legend(frameon=True, loc=loc_legend)

########## 3. the dm
i_row = 2

XRANGE = [-0.009, 0.009]
XLABEL = r"$\Delta m$"
YTICK = [1, 2, 3, 4, 5]
binvalue = np.array(YTICK)
YTICKLABELS = [r'$0.1< z_{\rm B} \leq 0.3$', 
                r'$0.3< z_{\rm B} \leq 0.5$', 
                r'$0.5< z_{\rm B} \leq 0.7$', 
                r'$0.7< z_{\rm B} \leq 0.9$', 
                r'$0.9< z_{\rm B} \leq 1.2$']
YLIM = [0.5, 5.5]

for i_col, inpath_list_tmp in enumerate(inpath_dm_list):

    ax = axs[i_row, i_col]

    # get values and plot
    Npoints = len(inpath_list_tmp)
    for i_val, inpath in enumerate(inpath_list_tmp):
        try:
            data = pd.read_csv(inpath)
        except FileNotFoundError:
            continue

        # the first row is for whole
        mvalue = (data.loc[1:, 'm1'].values + data.loc[1:, 'm2'].values) / 2.
        print('dm value', mvalue)
        # error
        merror = (data.loc[1:, 'm1_err'].values + data.loc[1:, 'm2_err'].values) / 2.
        print('dm error', merror)
        del data

        ax.errorbar(mvalue, binvalue + (i_val-Npoints/2.) * 0.1, xerr=merror,
                        color=COLORs_dm[i_val], marker='o', markersize=4, elinewidth=1.5, 
                        ls='none', label=LABELs_dm[i_val])

    for ibin in range(5):
        ax.axhline(y=1.5+ibin, color='black', ls='-', lw=1)

    ax.axvline(x=0.0, ls='dotted', label=None, color='k', linewidth=1.5)

    if i_col == 0:
        ax.set_yticks(ticks=YTICK, labels=YTICKLABELS)
    else:
        ax.set_yticks([])

    ax.tick_params(axis='y', length=0, width=0)

    ax.set_xlim(XRANGE[0], XRANGE[1])
    ax.set_ylim(YLIM[0], YLIM[1])
    ax.set_xlabel(XLABEL)

    # invert y-axis
    ax.invert_yaxis()

if outpath == 'show':
    plt.show()
    plt.close()
else:
    plt.savefig(outpath, dpi=300)
    plt.close()
    plt.switch_backend(backend_orig)
    print("plot saved as", outpath)    
