# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-02-09 11:14:12
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-09 17:18:34

### correlation functions of SKiLLS to test the linear shear bias calibration

import os
import sys
import numpy as np
import pandas as pd 

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, LogLocator, NullFormatter, NullLocator
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = False
mpl.rcParams['ytick.right'] = False

plt.rcParams['font.family'] = 'serif'
# font size
plt.rc('font', size=9)
# tex
plt.rcParams["text.usetex"] = True


# measured correlation function
inpath_list = ['/disks/shear10/ssli/K1000CS/skills_v07D7p1_Outputs/CF_A12_goldclasses_m283p283/xi_combined.csv',
                '/disks/shear10/ssli/K1000CS/skills_v07D7p1_Outputs/CF_A12_goldclasses_p283p283/xi_combined.csv',
                '/disks/shear10/ssli/K1000CS/skills_v07D7p1_Outputs/CF_A12_goldclasses_p283m283/xi_combined.csv',
                '/disks/shear10/ssli/K1000CS/skills_v07D7p1_Outputs/CF_A12_goldclasses_m283m283/xi_combined.csv']

# the input shear amplitude
xi_p_norm = (0.04)**2

# the shear bias
bias_file = '../B_Shear_Bias/results/m_skills_v07D7p1_LF_321_kidsPhotometry_A12_gold.csv'
tmp = pd.read_csv(bias_file)
m1_list = tmp.loc[1:, 'm1'].values
print('m1', m1_list)
m2_list = tmp.loc[1:, 'm2'].values
print('m2', m2_list)
m_list = (m1_list + m2_list)/2.
print('m', m_list)
error_budget = [0.019, 0.008, 0.007, 0.006, 0.006]

# number of redshift bins
nzbins = 5

# combine all measures
## weighted by the measurement uncertainties
for i_path, inpath in enumerate(inpath_list):
    cata = pd.read_csv(inpath) 
    if i_path == 0 :
        ito = cata['ito'].values
        jto = cata['jto'].values
        theta = cata['theta'].values / np.square(cata['xierr_p'].values)
        xip = cata['xi_p'].values / np.square(cata['xierr_p'].values)
        xip_err = 1./np.square(cata['xierr_p'].values)
        del cata
    else:
        theta += cata['theta'].values / np.square(cata['xierr_p'].values)
        xip += cata['xi_p'].values / np.square(cata['xierr_p'].values)
        xip_err += 1./np.square(cata['xierr_p'].values)
        del cata
## the final results
theta /= xip_err
xip /= xip_err
xip_err = np.sqrt(1./xip_err)
## save as dataframe
cata = pd.DataFrame({ 'ito': ito,
                        'jto': jto,
                        'theta': theta, 
                        'xi_p': xip,
                        'xierr_p': xip_err})

# compare to the input shear
cata.loc[:, 'xi_p'] /= xi_p_norm
cata.loc[:, 'xierr_p'] /= np.abs(xi_p_norm)

# compare to the estimated shear bias
for i in range(nzbins):
    for j in range(i, nzbins):
        mask_tmp = (cata['ito']==i) & (cata['jto']==j)
        ## scale to zero
        cata.loc[mask_tmp, 'xi_p'] -= (1+m_list[i])*(1+m_list[j])
CATAs = [cata]
del cata, mask_tmp

# plot
TITLE = None
outpath = 'show'
outpath = './plots/CF4m.pdf'

usetex = True

LABELs = None
COLORs = ['blue']
POINTs = ['o']
POINTSs = [2]

XRANGE = [0.5, 300]

timesX = False
yscaling = 100
YLABEL = r'$\Delta m_{\xi}~[10^{-2}]$'

which_half = 'xip'
hlines = [0.0]

xlog = True

XLABEL = r'$\theta$ (arcmin)'

YRANGE = [-6, 7]
# YRANGE = None

# get the error ranges
fill_between_xs_list=[]
fill_between_yLows_list=[]
fill_between_yHighs_list=[]
fill_between_COLORs_list=[]
fill_between_alphas_list=[]
for i in range(nzbins):
    error_i = error_budget[i]
    m_i = m_list[i]
    for j in range(i, nzbins):
        error_j = error_budget[j]
        m_j = m_list[j]

        error_range = (((1+m_i) * error_j)**2 + ((1+m_j) * error_i)**2)**0.5

        if XRANGE is not None:
            fill_between_xs_list.append([[XRANGE[0], XRANGE[1]]])
        else:
            fill_between_xs_list.append([[np.min(CATAs[0]['theta'].values), np.max(CATAs[0]['theta'].values)]])
        fill_between_yLows_list.append([-error_range])
        fill_between_yHighs_list.append([error_range])
        fill_between_COLORs_list.append(['gray'])
        fill_between_alphas_list.append([0.2])

########## the plot function
def _vhlines(vORh, lines, line_styles=None, line_colors=None, line_labels=None, line_widths=None, ax=plt):
    """
    Add vertical or horizontal lines to the main plots

    Parameters
    ----------
    vORh : {'v', 'h'}
        vertical line (v) or horizontal line (h)

    lines : array-like of floats
        where to place the lines (x-axis for vline, y-axis for hline)

    line_styles : array-like of linestyles, default: 'dashed'

    line_colors : array-like of colors, default: 'k'

    line_labels : array-like of strings, default: ''

    line_widths : array-like of floats, default: 1

    ax : matplotlib Axes object, default: matplotlib.pyplot

    Returns
    -------
        None
    """

    for i, line in enumerate(lines):

        if line_styles is not None:
            line_style = line_styles[i]
        else:
            line_style = '--'

        if line_colors is not None:
            line_color = line_colors[i]
        else:
            line_color = 'k'

        if line_labels is not None:
            line_label = line_labels[i]
        else:
            line_label = ''

        if line_widths is not None:
            line_width = line_widths[i]
        else:
            line_width = 1

        if vORh == 'v':
            ax.axvline(x=line, ls=line_style, label=line_label, color=line_color, linewidth=line_width)
        elif vORh == 'h':
            ax.axhline(y=line, ls=line_style, label=line_label, color=line_color, linewidth=line_width)
        else:
            raise Exception(f'Unsupported vORh value: {vORh}')

def CorrPlotFunc_half(outpath, nzbins,
                    CATAs, which_half,
                    COLORs, 
                    yscaling = 1e4, 
                    LABELs=None, 
                    LINEs=None, LINEWs=None, 
                    POINTs=None, POINTSs=None, ERRORSIZEs=None,
                    subLABEL_locX=0.1, subLABEL_locY=0.8,
                    XRANGE=None, YRANGE=None,
                    XLABEL=None, YLABEL=None, TITLE=None,
                    xtick_min_label=True, xtick_spe=None, ytick_min_label=True, ytick_spe=None,
                    vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                    hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                    xlog=False, invertX=False, ylog=False, invertY=False, 
                    loc_legend='best', legend_frame=False,
                    fill_between_xs_list=None, 
                    fill_between_yLows_list=None, fill_between_yHighs_list=None,
                    fill_between_COLORs_list=None, fill_between_alphas_list=None,
                    FIGSIZE=[6.4, 4.8],
                    TIGHT=False, 
                    N_rows=None, N_cols=None,
                    timesX=True,
                    xi_p_norm=1, xi_m_norm=1,
                    xi_p_err_col='xierr_p', xi_m_err_col='xierr_m'):
    """
    plot the correlation function
        only half: which_half = 'xip' or 'xim'
    used columns in cata:
        ito, jto, theta, xi_p, xi_m, xi_p_err_col, xi_m_err_col
    """

    if outpath != 'show':
        backend_orig = plt.get_backend()
        plt.switch_backend("agg")

    fig, ax = plt.subplots(nzbins, nzbins, sharex=True, sharey=True, figsize=FIGSIZE)
    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0)

    if TIGHT:
        # for x,y label
        fig.add_subplot(111, frameon=False)
        # hide tick and tick label of the big axes
        plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
        plt.grid(False)
        plt.xlabel(XLABEL)
        plt.ylabel(YLABEL)
        if TITLE is not None:
            plt.title(TITLE)

    # blank region
    for row_tmp in range(nzbins):
        for col_tmp in range(row_tmp+1, nzbins):
            ax[row_tmp, col_tmp].axis('off')

    # plot
    handles = []
    i_plot = 0
    for i in range(nzbins):
        for j in range(i, nzbins):
            row = j
            col = i

            for i_cata, cata in enumerate(CATAs):

                # get plot properties
                CR = COLORs[i_cata]
                if LABELs is not None:
                    LAB = LABELs[i_cata]
                else:
                    LAB = None
                if LINEs is not None:
                    LN = LINEs[i_cata]
                else:
                    LN = '--'
                if LINEWs is not None:
                    LW = LINEWs[i_cata]
                else:
                    LW = 1
                if POINTs is not None:
                    PI = POINTs[i_cata]
                else:
                    PI = 'o'
                if POINTSs is not None:
                    MS = POINTSs[i_cata]
                else:
                    MS = 2
                if ERRORSIZEs is not None:
                    ERRORSIZE = ERRORSIZEs[i_cata]
                else:
                    ERRORSIZE = 2

                # get values
                mask_tmp = (cata['ito']==i) & (cata['jto']==j)
                theta = cata.loc[mask_tmp, 'theta'].values
                if which_half == 'xip':
                    xi = cata.loc[mask_tmp, 'xi_p'].values / xi_p_norm
                    xierr = cata.loc[mask_tmp, xi_p_err_col].values / np.abs(xi_p_norm)
                elif which_half == 'xim':
                    xi = cata.loc[mask_tmp, 'xi_m'].values / xi_m_norm
                    xierr = cata.loc[mask_tmp, xi_m_err_col].values / np.abs(xi_m_norm)
                else:
                    raise Exception(f'Unsupported which_half {which_half}')
                del mask_tmp, cata

                # plot
                if timesX:
                    tmp = ax[row, col].errorbar(theta, theta*xi*yscaling, yerr=np.abs(theta*yscaling)*xierr,
                                                        color=CR, 
                                                        label=LAB, 
                                                        linestyle=LN, linewidth=LW, 
                                                        marker=PI, markersize=MS, capsize=ERRORSIZE)
                else:
                    tmp = ax[row, col].errorbar(theta, xi*yscaling, yerr=xierr*np.abs(yscaling),
                                                        color=CR, 
                                                        label=LAB, 
                                                        linestyle=LN, linewidth=LW, 
                                                        marker=PI, markersize=MS, capsize=ERRORSIZE)
                if (i_plot==0):
                    handles.append(tmp[0])
                del tmp

            if fill_between_xs_list is not None:
                fill_between_xs = fill_between_xs_list[i_plot]
                fill_between_yLows = fill_between_yLows_list[i_plot]
                fill_between_yHighs = fill_between_yHighs_list[i_plot]
                fill_between_COLORs = fill_between_COLORs_list[i_plot]
                fill_between_alphas = fill_between_alphas_list[i_plot]
                for i_fill, fill_between_x in enumerate(fill_between_xs):
                    ax[row, col].fill_between(fill_between_x, 
                                fill_between_yLows[i_fill]*yscaling, fill_between_yHighs[i_fill]*yscaling,
                                alpha=fill_between_alphas[i_fill],
                                color=fill_between_COLORs[i_fill])
            i_plot += 1

            LABEL = f'{i+1}-{j+1}'
            ax[row, col].text(subLABEL_locX, subLABEL_locY, LABEL, transform=ax[row, col].transAxes)

            if XRANGE is not None:
                ax[row, col].set_xlim(XRANGE[0], XRANGE[1])
            if YRANGE is not None:
                ax[row, col].set_ylim(YRANGE[0], YRANGE[1])

            if vlines is not None:
                _vhlines('v', vlines, line_styles=vline_styles, line_colors=vline_colors, line_labels=vline_labels, line_widths=vline_widths, ax=ax[row, col])
            if hlines is not None:
                _vhlines('h', hlines, line_styles=hline_styles, line_colors=hline_colors, line_labels=hline_labels, line_widths=hline_widths, ax=ax[row, col])

            if xlog:
                ax[row, col].set_xscale('log')
                ax[row, col].xaxis.set_minor_formatter(NullFormatter())
            if ylog:
                ax[row, col].set_yscale('log')
                ax[row, col].yaxis.set_minor_formatter(NullFormatter())

            if xtick_min_label:
                if xlog:
                    ax[row, col].xaxis.set_minor_locator(LogLocator(base=10.0, subs=None, numticks=10))
                else:
                    ax[row, col].xaxis.set_minor_locator(AutoMinorLocator())
            if ytick_min_label:
                if ylog:
                    ax[row, col].yaxis.set_minor_locator(LogLocator(base=10.0, subs=None, numticks=10))
                else:
                    ax[row, col].yaxis.set_minor_locator(AutoMinorLocator())

            if xtick_spe is not None:
                ax[row, col].set_xticks(xtick_spe[0])
                ax[row, col].set_xticklabels(xtick_spe[1])
            if ytick_spe is not None:
                ax[row, col].set_yticks(ytick_spe[0])
                ax[row, col].set_yticklabels(ytick_spe[1])

            if invertY:
                plt.gca().invert_yaxis()
            if invertX:
                plt.gca().invert_xaxis()

    if not TIGHT:
        fig.text(0.5, 0.03, XLABEL, ha='center')

        fig.text(0.05, 0.5, YLABEL, 
             horizontalalignment='center', verticalalignment='center',
             rotation='vertical')

        if TITLE is not None:
            fig.text(0.5, 0.90, TITLE, ha='center')

    if LABELs is not None:
        fig.legend(handles, LABELs, 
                loc = 'upper right',
                bbox_to_anchor=(0.92, 0.92), fancybox=False, shadow=False, frameon=False)

    if TIGHT:
        plt.tight_layout()

    if outpath == 'show':
        plt.show()
        plt.close()
    else:
        plt.savefig(outpath, dpi=300)
        plt.close()
        plt.switch_backend(backend_orig)
        print("Errorbar plot saved as", outpath)


# plot
CorrPlotFunc_half(outpath, nzbins,
                    CATAs, which_half,
                    COLORs, 
                    yscaling = yscaling, 
                    LABELs=LABELs, 
                    LINEs=None, LINEWs=None, 
                    POINTs=POINTs, POINTSs=POINTSs, ERRORSIZEs=None,
                    subLABEL_locX=0.65, subLABEL_locY=0.8,
                    XRANGE=XRANGE, YRANGE=YRANGE,
                    XLABEL=XLABEL, YLABEL=YLABEL, TITLE=TITLE,
                    xtick_min_label=True, xtick_spe=None, ytick_min_label=True, ytick_spe=None,
                    vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                    hlines=hlines, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                    xlog=xlog, invertX=False, ylog=False, invertY=False, 
                    loc_legend='best', legend_frame=False,
                    fill_between_xs_list=fill_between_xs_list, 
                    fill_between_yLows_list=fill_between_yLows_list, fill_between_yHighs_list=fill_between_yHighs_list,
                    fill_between_COLORs_list=fill_between_COLORs_list, fill_between_alphas_list=fill_between_alphas_list,
                    FIGSIZE=[6.4, 4.8],
                    TIGHT=False, 
                    N_rows=None, N_cols=None,
                    timesX=timesX)
