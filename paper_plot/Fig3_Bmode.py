# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-03-09 18:38:20
# @Last Modified by:   lshuns
# @Last Modified time: 2023-06-02 15:33:28

##### B mode in COSEBIs

import numpy as np
from scipy.stats import chi2

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rcParams
from mpl_toolkits.axes_grid1 import ImageGrid

##################### I/O

outpath = './plots/Bmode.pdf'
outpath = 'show'

# # ### 309c 
file_cosebis_single_309c = '/disks/shear10/ssli/Cat_to_Obs_K1000_P1-master/data/kids/\
cosebis_K1000_ALL_BLIND_C_bmodes_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid_nbins_theta_0.5_300.asc'
covdat_309c = '/disks/shear10/ssli/Cat_to_Obs_K1000_P1-master/data/covariance/outputs/\
Covariance_blindC_nMaximum_20_0.50_300.00_nBins5_NoiseOnly.ascii'

### 321
file_cosebis_321 = '/disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/COSEBIS/nmodes20/\
Bn_COSEBIS_K1000_ALL_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_theta_0.5_300_zbins_{:}_{:}.asc'
covdat_321 = '/disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/COV/output/\
Covariance_blindNONE_nMaximum_20_0.50_300.00_nBins5_NoiseOnly.ascii'

##################### workhorse

# Some font setting
plt.rc('font', size=16)
plt.rcParams["text.usetex"] = True
plt.rcParams['font.family'] = 'serif'

# set up the figure grid
tick_spacing_y = 0.5
tick_spacing_x = 5
minor_tick_spacing_x = 1

fig = plt.figure(figsize=(8, 6))
grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                 nrows_ncols=(5,3),
                 axes_pad=0.0,
                 aspect=False
                 )

# number of tomographic bins, COSEBIS modes and the value of the modes
ntomobin=5
# number of modes in data file
nmodestot=20
# number of modes to analyse
# nmodes=20 # full set
nmodes=5

# value of the modes
n=np.arange(1, nmodes+1, dtype=int)

# before we read in the per tomo bin combination data, we need to read in the full covariance from the mocks
## 309c
cov_309c = np.loadtxt(covdat_309c)
## 321
cov_321 = np.loadtxt(covdat_321)

# load Bn data
if 'file_cosebis_single_309c' in locals():
    Bn_in_309c = np.loadtxt(file_cosebis_single_309c)

# tomake different sub plots we count the grid square that we want to plot in
#initialising the counter
gridpos=-1

# read in B mode data per tomo bin combination
handles = []
LABELs = []
for iz in range(1,ntomobin+1):
    for jz in range(iz,ntomobin+1):

        gridpos=gridpos + 1
        ax=grid[gridpos]

        tomochar='%s_%s'%(iz,jz)

        ### 309c
        ipos = gridpos*nmodestot
        if 'file_cosebis_single_309c' in locals():
            Bn = Bn_in_309c[ipos:ipos+nmodes]
        else:
            Bnfile = file_cosebis_309c.format(iz, jz)
            Bn_in = np.loadtxt(Bnfile)
            Bn = Bn_in[0:nmodes]
        # breaking up the large covariance matrix to find the significance per tomographic bin combination
        cov_izjz = cov_309c[ipos:ipos+nmodes, ipos:ipos+nmodes]
        diagerr = np.sqrt(np.diagonal(cov_izjz))
        invcov = np.linalg.inv(cov_izjz)

        if gridpos==0:
            Bn_all_309c = Bn
        else:
            Bn_all_309c = np.append(Bn_prev_309c, Bn)
        Bn_prev_309c = Bn_all_309c

        # plot the results with annotations of the bin combination and p-value
        tmp = ax.errorbar(n, Bn*1e9, yerr=diagerr*1e9, 
            fmt='D', color='gray', markerfacecolor='none', markersize=4.0)
        if gridpos == 0:
            handles.append(tmp)
            LABELs.append('KiDS-1000-v1')

        # calculate the null chi-sq value and associated p-value
        chisq_null = np.matmul(Bn, (np.matmul(invcov, Bn)))
        pval_309c = chi2.sf(chisq_null, nmodes)

        ### 321
        Bnfile = file_cosebis_321.format(iz, jz)
        Bn_in = np.loadtxt(Bnfile)
        Bn = Bn_in[0:nmodes]
        # breaking up the large covariance matrix to find the significance per tomographic bin combination
        ipos = gridpos*nmodestot
        cov_izjz = cov_321[ipos:ipos+nmodes, ipos:ipos+nmodes]
        diagerr = np.sqrt(np.diagonal(cov_izjz))
        invcov = np.linalg.inv(cov_izjz)

        if gridpos==0:
            Bn_all_321 = Bn
        else:
            Bn_all_321 = np.append(Bn_prev_321, Bn)
        Bn_prev_321 = Bn_all_321

        # plot the results with annotations of the bin combination and p-value
        tmp = ax.errorbar(n, Bn*1e9, yerr=diagerr*1e9, 
            fmt='o', color='red', markerfacecolor='none', markersize=4.0)
        if gridpos == 0:
            handles.append(tmp)
            LABELs.append('KiDS-1000-v2')

        # calculate the null chi-sq value and associated p-value
        chisq_null = np.matmul(Bn, (np.matmul(invcov, Bn)))
        pval_321 = chi2.sf(chisq_null, nmodes)

        # add p values
        # pvalchar = 'p={:.2f} ({:.2f})'.format(pval_321, pval_309c)
        pvalchar = r'$p$' +' = {:.2f}'.format(pval_321)
        ax.annotate(pvalchar, xy=(0.95,0.9),xycoords='axes fraction',
            size=14, ha='right', va='top')

        # some general info
        ax.axhline(y=0, color='black', ls=':')
        ax.set_ylim(-1.7,2.1)
        ax.annotate(tomochar, xy=(0.2,0.9),xycoords='axes fraction',
            size=14, ha='right', va='top')
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(tick_spacing_y))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(minor_tick_spacing_x))
        ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing_x))

# add labels
fig.legend(handles, LABELs,
        loc = 'center', ncol=3,
        bbox_to_anchor=(0.5, 0.96), fancybox=True, shadow=True)

grid[6].set_ylabel(r'$B_n~[10^{-9}~{\rm rad}^2]$')
grid[12].set_xlabel(r'$n$')
grid[13].set_xlabel(r'$n$')
grid[14].set_xlabel(r'$n$')

# plt.tight_layout()
fig.subplots_adjust(top=0.91)

if outpath == 'show':
    plt.show()
    plt.close()
else:
    plt.savefig(outpath, dpi=300)
    plt.close()
    print("plot saved as", outpath)