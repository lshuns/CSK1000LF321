# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-05-16 18:14:00
# @Last Modified by:   lshuns
# @Last Modified time: 2023-09-11 13:48:55

### plot S8 constraints compared to other surveys 

import os

import numpy as np
import pandas as pd

from matplotlib.pyplot import cm
import matplotlib as mpl
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator, LogLocator, NullFormatter, NullLocator

# +++ general settings for plot
# mpl.use('Agg')
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
plt.rc('font', size=16)
plt.rcParams["text.usetex"] = True
plt.rcParams['font.family'] = 'serif'

############### I/O

## KiDS-1000-LF321 results
# MAP, Mean, max
# value, error low, error high
S8_us_LABEL = r'KiDS-1000-v2 (this work)'
S8_us_c = 'k'
# S8_us_list = [[0.774, 0.023, 0.027], 
#                 [0.765, 0.017, 0.028], 
#                 [0.771, 0.027, 0.023]] # Multinest Fiducial
S8_us_list = [[0.776, 0.027, 0.029], 
                [0.765, 0.023, 0.029], 
                [0.769, 0.029, 0.027]] # PolyChord Fiducial
S8_us_symbol = ['D', 's', 'o']

## results from other surveys
LABELs1 = ['DES Y3: Fiducial', 
            r'HSC Y3: $\xi$', 
            r'HSC Y3: $C_{\ell}$'
            ]
COLORs1 = ['blue',
        'chocolate', 'orange']
SYMBOLs1 = ['s', 
            'o', 'o']
S8_list_others1 = [[0.759, 0.023, 0.025],
                    [0.769, 0.034, 0.031],
                    [0.776, 0.033, 0.032]
                    ]

## results from joint analysis
LABELs2 = [
            'DES Y3+KiDS-1000-v1 Hybrid',
            'DES Y3 Hybrid',
            'KiDS-1000-v1 Hybrid'            
            ]
COLORs2 = ['m',
        'green', 'purple']
SYMBOLs2 = ['D', 'D', 'D']
S8_list_others2 = [[0.801, 0.023, 0.011],
                    [0.816, 0.028, 0.015],
                    [0.776, 0.027, 0.029]
                    ]

## results from Planck
LABEL_p = r'\textit{Planck}-2018'
COLOR_p = 'r'
SYMBOL_p = 's'
S8_list_p = [0.834, 0.016, 0.016]

# +++ plot related
outpath = 'show'
outpath = './plots/S8_KiDS_others.pdf'

XLIM = [0.6, 0.86]
MS = 7

YTICK = np.arange(len(COLORs1)+len(COLORs2)+2)
YTICKLABELS = [S8_us_LABEL] + LABELs1 + LABELs2 + [LABEL_p]
XLABEL = r"$S_8\equiv\sigma_8(\Omega_{\rm m}/0.3)^{0.5}$"
YLIM = [-0.5, len(YTICK) - 0.5]

# +++ plot
fig, ax = plt.subplots()

# plot our results
for i_val, val_us in enumerate(S8_us_list):

    yval = (i_val - 1) * 0.2
    plt.errorbar(val_us[0], yval, xerr=[[val_us[1]], [val_us[2]]],
                    color=S8_us_c, marker=S8_us_symbol[i_val], markersize=MS, ls='none')

    if i_val == 0:
        plt.fill_between(np.array([val_us[0] - val_us[1], 
                        val_us[0] + val_us[2]]), 
                    YLIM[0], YLIM[1], 
                    edgecolor=S8_us_c, facecolor=S8_us_c, alpha=0.1)

# plot other results
for i_val, val in enumerate(S8_list_others1):

    yval = i_val + 1

    CR = COLORs1[i_val]
    MK = SYMBOLs1[i_val]

    plt.errorbar(val[0], yval, xerr=[[val[1]], [val[2]]],
                    color=CR, marker=MK, markersize=MS, ls='none')
## a dotted line
plt.hlines(yval+0.5, XLIM[0], XLIM[1], colors='gray', linestyles='dotted')

# plot hybrid results
for i_val, val in enumerate(S8_list_others2):

    yval2 = yval + i_val + 1

    CR = COLORs2[i_val]
    MK = SYMBOLs2[i_val]

    plt.errorbar(val[0], yval2, xerr=[[val[1]], [val[2]]],
                    color=CR, marker=MK, markersize=MS, ls='none')
## a dotted line
plt.hlines(yval2+0.5, XLIM[0], XLIM[1], colors='gray', linestyles='dotted')

# plot Planck
plt.errorbar(S8_list_p[0], yval2+1, xerr=[[S8_list_p[1]], [S8_list_p[2]]],
                color=COLOR_p, marker=SYMBOL_p, markersize=MS, ls='none')

plt.xticks(ticks=[0.7, 0.75, 0.8, 0.85])
ax.xaxis.set_minor_locator(AutoMinorLocator())

plt.yticks(ticks=YTICK, labels=YTICKLABELS, horizontalalignment='left')
plt.tick_params(axis='y', length=0, width=0)

# Customize the yticklabels to be inside the plot
ax.tick_params(axis='y', pad=-20)

# Change ytick label colours
for color, label in zip([S8_us_c] + COLORs1 + COLORs2 + [COLOR_p], ax.get_yticklabels()):
    label.set_color(color)

plt.xlim(XLIM[0], XLIM[1])
plt.ylim(YLIM[0], YLIM[1])

xlabel = ax.set_xlabel(XLABEL)
xlabel.set_position((0.68, -0.32)) 

# invert y-axis
plt.gca().invert_yaxis()

plt.tight_layout()

if outpath == 'show':
    plt.show()
else:
    plt.savefig(outpath, dpi=300)
    print('plot saved in', outpath)
plt.close()

