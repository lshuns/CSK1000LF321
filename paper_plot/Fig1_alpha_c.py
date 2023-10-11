# -*- coding: utf-8 -*-
# @Author: ssli
# @Date:   2022-09-22 16:45:30
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-31 17:13:00

### alpha and c for tomographic bins

import pandas as pd 
import numpy as np

import plotting

# >>>>>>> I/O

inpath_list = ['../B_Shear_Bias/results/alpha_K1000_LF309c_public.csv',
                '../B_Shear_Bias/results/alpha_K1000_LF321_alphaRecal.csv']
COLORs = ['gray', 'red']
LABELs = ['KiDS-1000-v1', 'KiDS-1000-v2']
LINEs = ['', '']
POINTs = ['x', 'o']
POINTSs = [4, 6]

# >>>>>>> general info
outpath = 'show'
outpath = './plots/alpha_c.pdf'

usetex = True

N_plots = 4
N_rows = 2
N_cols = 2
sharex = True
sharey = False

subLABEL_bbox = dict(facecolor='white', edgecolor='black', boxstyle='round')

XLABEL = r'$z_{\rm B}$'
XRANGE = [0.1, 1.2]

yscaling_list = [1e2, 1e4]
YLABEL_list = [r'$\alpha~[10^{-2}]$', None, r'$c~[10^{-4}]$', None]
YRANGE_list = [[-6, 6], [-6, 6], [-9, 9], [-9, 9]]
no_xticklabels_list = None
no_yticklabels_list = [False, True, False, True]
subLABEL_list = [r'$\epsilon_1$', r'$\epsilon_2$', r'$\epsilon_1$', r'$\epsilon_2$']

LABEL_position = 'top'
LABEL_cols = 2

font_size = 16

hlines = [0.0]
hline_styles = ['--']
hline_colors = ['k']

# >>>>>>> plot

# get values
cata_list = []
for inpath in inpath_list:
    cata = pd.read_csv(inpath)
    cata_list.append(cata)
    del cata

# collect values for plot
xvals_list = []
yvals_list = []
yerrs_list = []
fill_between_xs_list=[] 
fill_between_yLows_list=[] 
fill_between_yHighs_list=[]
fill_between_COLORs_list=[]
fill_between_alphas_list=[]
COLORs_list = []
LABELs_list = []
LINEs_list = []
POINTs_list = []
POINTSs_list = []
for i_col, col_tmp in enumerate(['alpha', 'c']):

    # loop over component
    for e_col in ['1', '2']:

        # collect values from each cata
        xvals = []
        yvals = []
        yerrs = []
        fill_between_xs=[] 
        fill_between_yLows=[] 
        fill_between_yHighs=[]
        fill_between_COLORs=COLORs
        fill_between_alphas=[0.3, 0.3, 0.3, 0.3]
        for i_cata, cata in enumerate(cata_list):

            ## the first is for whole
            fill_between_xs.append([np.min(cata.loc[1:, 'Z_B_min'].values), np.max(cata.loc[1:, 'Z_B_max'].values)])
            fill_between_yLows.append(cata[f'{col_tmp}{e_col}'].values[0]* yscaling_list[i_col] - cata[f'{col_tmp}{e_col}_err'].values[0]* yscaling_list[i_col])
            fill_between_yHighs.append(cata[f'{col_tmp}{e_col}'].values[0]* yscaling_list[i_col] + cata[f'{col_tmp}{e_col}_err'].values[0]* yscaling_list[i_col])

            cata = cata.loc[1:, :]
            # center ZB as x
            xvals.append((cata['Z_B_min'].values + cata['Z_B_max'].values)/2.-0.005*(i_cata+1))
            # y from bias values
            yvals.append(cata[f'{col_tmp}{e_col}'].values* yscaling_list[i_col])

            yerrs.append([cata[f'{col_tmp}{e_col}_err'].values* yscaling_list[i_col], cata[f'{col_tmp}{e_col}_err'].values* yscaling_list[i_col]])

        xvals_list.append(xvals)
        yvals_list.append(yvals)
        yerrs_list.append(yerrs)
        fill_between_xs_list.append(fill_between_xs)
        fill_between_yLows_list.append(fill_between_yLows)
        fill_between_yHighs_list.append(fill_between_yHighs)
        fill_between_COLORs_list.append(fill_between_COLORs)
        fill_between_alphas_list.append(fill_between_alphas)
        COLORs_list.append(COLORs)
        LABELs_list.append(LABELs)
        LINEs_list.append(LINEs)
        POINTs_list.append(POINTs)
        POINTSs_list.append(POINTSs)

plotting.ErrorPlotFunc_subplots(outpath, N_plots,
                            xvals_list, yvals_list, yerrs_list,
                            COLORs_list, LABELs_list=LABELs_list, LINEs_list=LINEs_list, LINEWs_list=None, 
                            POINTs_list=POINTs_list, POINTSs_list=POINTSs_list, ERRORSIZEs_list=None,
                            subLABEL_list=subLABEL_list, subLABEL_locX=0.1, subLABEL_locY=0.8, subLABEL_bbox=subLABEL_bbox,
                            XRANGE=XRANGE, YRANGE=None,
                            XLABEL=XLABEL, YLABEL=None, TITLE=None,
                            xtick_min_label=True, xtick_spe=None, ytick_min_label=True, ytick_spe=None,
                            vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                            hlines=hlines, hline_styles=hline_styles, hline_colors=hline_colors, hline_labels=None, hline_widths=None,
                            xlog=False, invertX=False, ylog=False, invertY=False, 
                            loc_legend='best', legend_frame=True,
                            font_size=font_size, usetex=usetex,
                            fill_between_xs_list=fill_between_xs_list, 
                            fill_between_yLows_list=fill_between_yLows_list, fill_between_yHighs_list=fill_between_yHighs_list,
                            fill_between_COLORs_list=fill_between_COLORs_list, fill_between_alphas_list=fill_between_alphas_list,
                            LABEL_position=LABEL_position, LABEL_position_SUBid=0,
                            LABEL_cols=LABEL_cols,
                            FIGSIZE=[6.4, 4.8],
                            TIGHT=False, 
                            N_rows=N_rows, N_cols=N_cols,
                            sharex=sharex, sharey=sharey,
                            YLABEL_list = YLABEL_list,
                            YRANGE_list = YRANGE_list,
                            no_xticklabels_list=no_xticklabels_list,
                            no_yticklabels_list=no_yticklabels_list)