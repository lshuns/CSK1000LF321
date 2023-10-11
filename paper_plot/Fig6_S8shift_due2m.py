# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-05-09 16:00:49
# @Last Modified by:   lshuns
# @Last Modified time: 2023-06-01 09:24:01

### plot S8 shift due to the different m scenarios

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

# ++++++++++++++++++++++ I/O
# run_name = 'multinest_nIA'
# err_range = [-0.023, 0.027]

run_name = 'polychord_nIA'
err_range = [-0.027, 0.029]
err_colour = 'gray'

m_variations = ['no_m_bias', 'with_raw_m_bias',
                'test_m_sizeU', 'test_m_sizeD',
                'test_m_qU', 'test_m_qD',
                'test_m_nU', 'test_m_nD']
LABELs = [r'No $m$ correction', r'Using $m_{\rm raw}$',
        r'Input size shift up', r'Input size shift down',
        r'Input $q$ shift up', r'Input $q$ shift down', 
        r'Input $n$ shift up', r'Input $n$ shift down']
COLORs = ['slategrey', 'darkslategray',
        'orange', 'blue',
        'm', 'c',
        'purple', 'green']
SYMBOLs = ['.', 'o', 
            '>', '<',
            '>', '<',
            '>', '<']

indir = f'/disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/{run_name}'
# results with the fiducial m
infile0 = f'maxpost_{run_name}_start.txt'
# results for comparison
infile_list = [f'maxpost_{run_name}_{m_variation}.txt' for m_variation in m_variations]

S8_col_id = 12

# +++ sigma
sigma = (err_range[1]-err_range[0])/2.
print(">>>> sigma", sigma)

# +++ plot related
outpath = 'show'
outpath = './plots/S8shift_due2m.pdf'

Nval = len(COLORs)
XLIM = [-0.04, 0.02]
MSs = [7] * Nval

YTICK = np.arange(Nval)
YTICKLABELS = LABELs
XLABEL = r"$\Delta S_8$"
YLIM = [-0.5, Nval - 0.5]

# +++ plot
fig, ax = plt.subplots()

# load the fiducial results
S80 = np.loadtxt(os.path.join(indir, infile0))[-1, S8_col_id]
print(">>>>> check if S8 make sense", S80)

xmin = 0
xmax = 0
for i, infile in enumerate(infile_list):

    S8 = np.loadtxt(os.path.join(indir, infile))[-1, S8_col_id]

    xval = S8 - S80
    print("+++ shift in sigma", infile, xval / sigma)

    if i > 0:
        if xval > xmax:
            xmax = xval
        if xval < xmin:
            xmin = xval

    yval = i

    CR = COLORs[i]
    MK = SYMBOLs[i]
    MS = MSs[i]

    plt.errorbar([xval], [yval],
                    color=CR, marker=MK, markersize=MS, ls='none')

plt.axvline(x=0, color='k', ls = '-', lw=1)

#### the sigma region
err_ratios = [0.5, 0.3, 0.1]
alphas = [0.1, 0.3, 0.5]
for i_ratio, err_ratio in enumerate(err_ratios):
    plt.fill_between(np.array([err_range[0]*err_ratio, 
                            err_range[1]*err_ratio]), 
                        YLIM[0], YLIM[1], 
                        edgecolor=err_colour, facecolor=err_colour, alpha=alphas[i_ratio])

plt.axvline(x=xmin, color='k', ls = '--', lw=1)
plt.axvline(x=xmax, color='k', ls = '--', lw=1)
print(">>>>> S8 shift", xmin, xmax)

plt.xticks(ticks=[-0.02, 0.00, 0.02])
ax.xaxis.set_minor_locator(AutoMinorLocator())

plt.yticks(ticks=YTICK, labels=YTICKLABELS, horizontalalignment='left')
plt.tick_params(axis='y', length=0, width=0)

# Customize the yticklabels to be inside the plot
ax.tick_params(axis='y', pad=-15)

# Change ytick label colours
for color, label in zip(COLORs, ax.get_yticklabels()):
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

############### multinest
# >>>>> check if S8 make sense 0.77367
# >>>>> S8 shift -0.0025999999999999357 0.0022400000000000198
# plot saved in ./plots/S8shift_due2m.pdf
# Elapsed:0:06.46,User=6.736,System=5.201,CPU=184.6%.

# ############### polychord
# >>>> sigma 0.028
# >>>>> check if S8 make sense 0.77564
# +++ shift in sigma maxpost_polychord_nIA_no_m_bias.txt 0.614285714285714
# +++ shift in sigma maxpost_polychord_nIA_with_raw_m_bias.txt -0.03107142857142991
# +++ shift in sigma maxpost_polychord_nIA_test_m_sizeU.txt -0.07857142857142785
# +++ shift in sigma maxpost_polychord_nIA_test_m_sizeD.txt 0.05821428571428644
# +++ shift in sigma maxpost_polychord_nIA_test_m_qU.txt -0.10499999999999796
# +++ shift in sigma maxpost_polychord_nIA_test_m_qD.txt 0.05857142857142767
# +++ shift in sigma maxpost_polychord_nIA_test_m_nU.txt 0.03749999999999983
# +++ shift in sigma maxpost_polychord_nIA_test_m_nD.txt -0.05464285714285826
# >>>>> S8 shift -0.0029399999999999427 0.0016399999999999748
# plot saved in ./plots/S8shift_due2m.pdf
