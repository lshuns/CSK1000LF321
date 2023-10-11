# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-05-19 16:56:54
# @Last Modified by:   lshuns
# @Last Modified time: 2023-06-01 09:19:17

### plot S8 and nuisance para contours with chaincosumer

import numpy as np
from chainconsumer import ChainConsumer

from matplotlib.pyplot import cm
import matplotlib as mpl
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator, LogLocator, NullFormatter, NullLocator

# +++ general settings for plot
# mpl.use('Agg')
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = False
mpl.rcParams['ytick.right'] = False
plt.rc('font', size=12)
plt.rcParams["text.usetex"] = True
plt.rcParams['font.family'] = 'serif'

# >>>>>>>>>>> I/O

outpath = './plots/S8_nuisance.pdf'
# outpath = 'show'

inpath_list = ['/disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/polychord_nIA/chain_polychord_nIA.txt',
                '/disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/multinest_nIA/output_multinest_nIA.txt']
LABELs = [r'KiDS-1000-v2: \textsc{PolyChord}',
            r'KiDS-1000-v2: \textsc{MultiNest}']
COLORs = ['k', 'grey']
ALPHAs = [0, 0]
LSs = ['-', '--']

paras = ['s8', 'a_ia', 
        'deltaz_1', 'deltaz_2', 'deltaz_3', 'deltaz_4', 'deltaz_5']
paras_label = [r'$S_8$', r'$A_{\rm IA}$',
                r'$\delta_{z,1}$', r'$\delta_{z,2}$', r'$\delta_{z,3}$', r'$\delta_{z,4}$', r'$\delta_{z,5}$'
            ]

paras_truth = [None, None, 0, -0.002, -0.013, -0.011, 0.006]

# >>>>>>>>>>>>>>>>>>>>> the parameter names in MCMC chain
input_names={ 'omch2' :'cosmological_parameters--omch2', 
              'ombh2':'cosmological_parameters--ombh2', 
              'h'    :'cosmological_parameters--h0', 
              'h_out'     :'COSMOLOGICAL_PARAMETERS--H0',
              'n_s'   :'cosmological_parameters--n_s', 
              's8_in':'cosmological_parameters--s_8_input', 
              'logt_agn':'halo_model_parameters--logt_agn',
              'a_bar'      :'halo_model_parameters--a', 
              'a_ia':'intrinsic_alignment_parameters--a', 
              'alpha_ia': 'intrinsic_alignment_parameters--alpha',
              'deltaz_uncorr_1' :'nofz_shifts--uncorr_bias_1', 
              'deltaz_uncorr_2' :'nofz_shifts--uncorr_bias_2', 
              'deltaz_uncorr_3' :'nofz_shifts--uncorr_bias_3', 
              'deltaz_uncorr_4' :'nofz_shifts--uncorr_bias_4', 
              'deltaz_uncorr_5' :'nofz_shifts--uncorr_bias_5', 
              'uncorr_m1' :'shear_calibration_parameters--uncorr_m1', 
              'uncorr_m2' :'shear_calibration_parameters--uncorr_m2', 
              'uncorr_m3' :'shear_calibration_parameters--uncorr_m3', 
              'uncorr_m4' :'shear_calibration_parameters--uncorr_m4', 
              'uncorr_m5' :'shear_calibration_parameters--uncorr_m5', 
              'delta_c' : 'shear_c_bias--delta_c',
              's8':    'COSMOLOGICAL_PARAMETERS--S_8', 
              'sigma_8':  'COSMOLOGICAL_PARAMETERS--SIGMA_8', 
              'a_s': 'COSMOLOGICAL_PARAMETERS--A_S', 
              'omega_m':    'COSMOLOGICAL_PARAMETERS--OMEGA_M', 
              'omega_nu':    'COSMOLOGICAL_PARAMETERS--OMEGA_NU', 
              'omega_lambda':       'COSMOLOGICAL_PARAMETERS--OMEGA_LAMBDA', 
              'theta_mc':      'COSMOLOGICAL_PARAMETERS--COSMOMC_THETA', 
              'deltaz_1':'NOFZ_SHIFTS--BIAS_1', 
              'deltaz_2':'NOFZ_SHIFTS--BIAS_2', 
              'deltaz_3':'NOFZ_SHIFTS--BIAS_3', 
              'deltaz_4':'NOFZ_SHIFTS--BIAS_4', 
              'deltaz_5':'NOFZ_SHIFTS--BIAS_5', 
              'deltaz_out_1':'DELTA_Z_OUT--BIN_1', 
              'deltaz_out_2':'DELTA_Z_OUT--BIN_2', 
              'deltaz_out_3':'DELTA_Z_OUT--BIN_3', 
              'deltaz_out_4':'DELTA_Z_OUT--BIN_4', 
              'deltaz_out_5':'DELTA_Z_OUT--BIN_5',
              'm1':'SHEAR_CALIBRATION_PARAMETERS--M1',
              'm2':'SHEAR_CALIBRATION_PARAMETERS--M2',
              'm3':'SHEAR_CALIBRATION_PARAMETERS--M3',
              'm4':'SHEAR_CALIBRATION_PARAMETERS--M4',
              'm5':'SHEAR_CALIBRATION_PARAMETERS--M5',
              'prior':'prior', 
              'like':'like', 
              'post':'post', 
              'weight':'weight'}

# >>>>>>>>>>> workhorse

# load chains
chain_list = []
weight_list = []
for i_file, inpath in enumerate(inpath_list):

    # get column names
    with open(inpath) as tmp:
        cols_names = (tmp.readline().replace('#', '').lower()).split()
    cols_index = [cols_names.index(input_names[para].lower()) for para in paras]
    wei_index = cols_names.index('weight')
    
    # get values
    cata = np.loadtxt(inpath)
    chain = cata[:, cols_index]
    weight = cata[:, wei_index]
    del cata

    chain_list.append(chain)
    weight_list.append(weight)
    del chain

# build chains
c = ChainConsumer()
for i_chain, chain in enumerate(chain_list):
    c.add_chain(chain, weights=weight_list[i_chain],
        linestyle=LSs[i_chain], color=COLORs[i_chain], 
        parameters=paras_label, name=LABELs[i_chain], shade_alpha=ALPHAs[i_chain])

# plot
c.configure(kde=1.5, plot_hists=True, shade_gradient=0.,
        diagonal_tick_labels=False,
        label_font_size=12, tick_font_size=12, 
        serif=True, legend_color_text=True, linewidths=1.5, statistics="max", shade=True,
        bar_shade=False,
        usetex=True)
c.configure_truth(linestyle=':')
fig = c.plotter.plot(truth=paras_truth)

# adjust each ax
axs = fig.get_axes()

xlims = [[0.65, 0.85], None, None, None, None, None, None,
        [0.65, 0.85], [-0.2, 1.1], None, None, None, None, None, 
        [0.65, 0.85], [-0.2, 1.1], [-0.04, 0.04], None, None, None, None, 
        [0.65, 0.85], [-0.2, 1.1], [-0.04, 0.04], [-0.04, 0.04], None, None, None,
        [0.65, 0.85], [-0.2, 1.1], [-0.04, 0.04], [-0.04, 0.04], [-0.04, 0.04], None, None,
        [0.65, 0.85], [-0.2, 1.1], [-0.04, 0.04], [-0.04, 0.04], [-0.04, 0.04], [-0.04, 0.04], None,
        [0.65, 0.85], [-0.2, 1.1], [-0.04, 0.04], [-0.04, 0.04], [-0.04, 0.04], [-0.04, 0.04], [-0.04, 0.04]
        ]

xticks = [[0.7, 0.8], None, None, None, None, None, None,
        [0.7, 0.8], [0, 1], None, None, None, None, None,
        [0.7, 0.8], [0, 1], [-0.02, 0.02], None, None, None, None,
        [0.7, 0.8], [0, 1], [-0.02, 0.02], [-0.02, 0.02], None, None, None,
        [0.7, 0.8], [0, 1], [-0.02, 0.02], [-0.02, 0.02], [-0.02, 0.02], None, None,
        [0.7, 0.8], [0, 1], [-0.02, 0.02], [-0.02, 0.02], [-0.02, 0.02], [-0.02, 0.02], None,
        [0.7, 0.8], [0, 1], [-0.02, 0.02], [-0.02, 0.02], [-0.02, 0.02], [-0.02, 0.02], [-0.02, 0.02]
        ]
xtick_labels = [['']*2, None, None, None, None, None, None,
        ['']*2, ['']*2, None, None, None, None, None,
        ['']*2, ['']*2, ['']*2, None, None, None, None,
        ['']*2, ['']*2, ['']*2, ['']*2, None, None, None,
        ['']*2, ['']*2, ['']*2, ['']*2, ['']*2, None, None,
        ['']*2, ['']*2, ['']*2, ['']*2, ['']*2, ['']*2, None,
        ['0.7', '0.8'], ['0', '1'], ['-0.02', '0.02'], ['-0.02', '0.02'], ['-0.02', '0.02'], ['-0.02', '0.02'], ['-0.02', '0.02']
        ]

ylims = [None, None, None, None, None, None, None,
        [-0.2, 1.1], None, None, None, None, None, None,
        [-0.04, 0.04], [-0.04, 0.04], None, None, None, None, None, 
        [-0.04, 0.04], [-0.04, 0.04], [-0.04, 0.04], None, None, None, None,
        [-0.04, 0.04], [-0.04, 0.04], [-0.04, 0.04], [-0.04, 0.04], None, None, None,
        [-0.04, 0.04], [-0.04, 0.04], [-0.04, 0.04], [-0.04, 0.04], [-0.04, 0.04], None, None,
        [-0.04, 0.04], [-0.04, 0.04], [-0.04, 0.04], [-0.04, 0.04], [-0.04, 0.04], [-0.04, 0.04], None
        ]
yticks = [None, None, None, None, None, None, None,
        [0, 1.], None, None, None, None, None, None,
        [-0.02, 0.02], [-0.02, 0.02], None, None, None, None, None, 
        [-0.02, 0.02], [-0.02, 0.02], [-0.02, 0.02], None, None, None, None,
        [-0.02, 0.02], [-0.02, 0.02], [-0.02, 0.02], [-0.02, 0.02], None, None, None,
        [-0.02, 0.02], [-0.02, 0.02], [-0.02, 0.02], [-0.02, 0.02], [-0.02, 0.02], None, None,
        [-0.02, 0.02], [-0.02, 0.02], [-0.02, 0.02], [-0.02, 0.02], [-0.02, 0.02], [-0.02, 0.02], None
        ]
ytick_labels = [None, None, None, None, None, None, None,
        ['0', '1'], None, None, None, None, None, None,
        ['-0.02', '0.02'], ['']*2, None, None, None, None, None, 
        ['-0.02', '0.02'], ['']*2, ['']*2, None, None, None, None,
        ['-0.02', '0.02'], ['']*2, ['']*2, ['']*2, None, None, None,
        ['-0.02', '0.02'], ['']*2, ['']*2, ['']*2, ['']*2, None, None,
        ['-0.02', '0.02'], ['']*2, ['']*2, ['']*2, ['']*2, ['']*2, None
        ]

for i_ax, ax in enumerate(axs):
    if xlims[i_ax] is not None:
        ax.set_xlim(xlims[i_ax][0], xlims[i_ax][1])
    if ylims[i_ax] is not None:
        ax.set_ylim(ylims[i_ax][0], ylims[i_ax][1])
    if xticks[i_ax] is not None:
        ax.set_xticks(ticks=xticks[i_ax], labels=xtick_labels[i_ax])
        ax.xaxis.set_minor_locator(AutoMinorLocator())
    if yticks[i_ax] is not None:
        ax.set_yticks(ticks=yticks[i_ax], labels=ytick_labels[i_ax])
        ax.yaxis.set_minor_locator(AutoMinorLocator())
    if (i_ax == 7) or (i_ax==14) or (i_ax==21) or (i_ax==28) or (i_ax==35) or (i_ax==42) or (i_ax==49):
        ylabel = ax.get_ylabel()
        ax.set_ylabel('')
        ax.text(-0.5, 0.5, ylabel, va='center', rotation='vertical', transform=ax.transAxes)

if outpath == 'show':
    plt.show()
    plt.close()
else:
    plt.savefig(outpath, dpi=300)
    plt.close()
    print("plot saved as", outpath)    