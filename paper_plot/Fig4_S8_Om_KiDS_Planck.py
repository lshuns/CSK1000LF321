# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-05-09 18:38:03
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-31 21:42:10

### plot S8 vs Omega_m contours with chaincosumer

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

outpath = './plots/S8_Om_KiDS_Planck.pdf'
# outpath = 'show'

inpath_list = ['/disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/polychord_nIA/chain_polychord_nIA.txt',
                '/disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/multinest_nIA/output_multinest_nIA.txt',
                '/disks/shear10/ssli/K1000CS/other_chains/planck2018/base_plikHM_TTTEEE_lowl_lowE.txt']
Om_col_id_list = [15, 15, 31]
S8_col_id_list = [12, 12, 36]
weight_col_id_list = [-1, -1, 0]
LABELs = [r'KiDS-1000-v2: \textsc{PolyChord}',
            r'KiDS-1000-v2: \textsc{MultiNest}',
            r'\textit{Planck}-2018']
COLORs = ['k', 'grey', 'r']
ALPHAs = [0, 0, 0]
LSs = ['-', '--', '-']

# >>>>>>>>>>> workhorse

# load chains
chain_list = []
weight_list = []
for i_file, inpath in enumerate(inpath_list):
    cata = np.loadtxt(inpath)

    chain = cata[:, [Om_col_id_list[i_file], S8_col_id_list[i_file]]]
    if weight_col_id_list[i_file] is not None:
        weight = cata[:, weight_col_id_list[i_file]]
    else:
        weight = None
    del cata

    chain_list.append(chain)
    weight_list.append(weight)
    del chain

# build chains
c = ChainConsumer()
for i_chain, chain in enumerate(chain_list):
    c.add_chain(chain, weights=weight_list[i_chain],
        linestyle=LSs[i_chain], color=COLORs[i_chain], 
        parameters=[r'$\Omega_m$', r'$S_8$'], name=LABELs[i_chain], shade_alpha=ALPHAs[i_chain])

# plot
c.configure(kde=1.5, plot_hists=False, shade_gradient=0.,
        diagonal_tick_labels=False,
        label_font_size=12, tick_font_size=12, 
        serif=True, legend_color_text=True, linewidths=1.5, statistics="max", shade=True,
        usetex=True)
fig = c.plotter.plot(figsize='column')

# change legend position
ax = fig.get_axes()[0]
ax.get_legend().remove()
for i_chain, chain_label in enumerate(LABELs):
    ax.text(0.95, 0.95 - i_chain * 0.05, chain_label, color=COLORs[i_chain], transform=ax.transAxes, ha='right')

plt.xlim(0.1, 0.6)
plt.ylim(0.68, 0.9)

plt.tight_layout()

if outpath == 'show':
    plt.show()
    plt.close()
else:
    plt.savefig(outpath, dpi=300)
    plt.close()
    print("plot saved as", outpath)    