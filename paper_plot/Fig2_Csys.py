# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-03-08 14:19:29
# @Last Modified by:   lshuns
# @Last Modified time: 2023-06-07 09:56:10

### compare the systematics to signal in two-point statistics

import os
import numpy as np
from astropy.io import fits

import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
from matplotlib import rcParams

# >>>>>>>>>>>>> I/O

# systematics
## 309c results
filetop_309c = '/disks/shear10/ssli/K1000CS/K1000_info/CSys/CSys_BLIND_C_5Z'
LFVER_309c = 'LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid'
## 321 results
filetop_321 = '/disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/Csys/CSys_BLIND_NONE_5Z'
LFVER_321 = 'LF_glab_321_v2_A12_goldclasses'

# folder with theory signal
folder_theory = '/disks/shear10/ssli/K1000CS/K1000_info/ForBG2/new_outputs/test_output_S8_fid_test/chain/output_test_A/shear_xi_plus_binned'

# file with covariance
# inpath_cov = '/disks/shear10/ssli/K1000CS/K1000_info/KiDS1000_cosmis_shear_data_release/data_fits/xipm_KIDS1000_BlindC_with_m_bias_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid.fits'
inpath_cov = '/disks/shear10/ssli/K1000CS/K1000_info/covariance_matrix_lf321.mat'

# output
outpath = 'show'
outpath = './plots/Csys.pdf'

# >>>>>>>>>>>>> general setups

plt.rc('font', size=16)
plt.rcParams["text.usetex"] = True
plt.rcParams['font.family'] = 'serif'

# set up the figure grid
tick_spacing = 2
fig,axes = plt.subplots(5,3,figsize=(8, 12),gridspec_kw={'hspace': 0, 'wspace': 0})

# number of tomographic bins
ntomobin = 5

# number of angular bins
nbins = 9

# initialising the counter
gridpos = -1

# >>>>>>>>>>>>> work horse

# read in the covariance
if '.fits' in os.path.basename(inpath_cov):
    with fits.open(inpath_cov) as f:
        COVMAT = f[1].data
elif ".mat" in os.path.basename(inpath_cov):
    COVMAT = np.loadtxt(inpath_cov)

# read in Csys data per tomo bin                                                                                                                                               
handles = []
LABELs = []
for iz in range(0,ntomobin):
    for jz in range(iz,ntomobin):

        gridpos = gridpos + 1
        
        # which grid cell do we want to plot this in?
        grid_x_E=int(gridpos/3)
        grid_y_E=gridpos % 3    #remainder
        ax=axes[grid_x_E,grid_y_E]

        # read in the expectation value for the cosmic shear signal
        xiptheory = np.loadtxt(os.path.join(folder_theory, 'bin_%d_%d.txt'%(jz+1,iz+1)))
        thetatheory = np.loadtxt(os.path.join(folder_theory, 'theta_bin_%d_%d.txt'%(jz+1,iz+1)))

        # we want to include a band which shows the error
        if COVMAT.shape[0] == 270:
            istart=gridpos  # because of only containing xi+ 
            diag_for_izjz=np.diag(COVMAT[istart:istart+nbins,istart:istart+nbins])
        elif COVMAT.shape[0] == 135:
            istart=gridpos*nbins  #factor of 2 because of xi+ and xi-
            diag_for_izjz=np.diag(COVMAT[istart:istart+nbins,istart:istart+nbins])
        fact = 0.1
        tmp = ax.fill_between(thetatheory, fact*np.sqrt(diag_for_izjz)/xiptheory *-1.0, fact*np.sqrt(diag_for_izjz)/xiptheory, 
                color='lightblue', linestyle=':', alpha=0.3)
        if gridpos == 0:
            handles.append(tmp)
            LABELs.append(r'$0.1\sigma_{\xi_+}$')

        # the dashed lines show 2 per cent 
        ax.axhline(y=0, color='black')
        ax.axhline(y=0.02, color='black', ls=':')
        ax.axhline(y=-0.02, color='black', ls=':')

        # read in the systematics
        tomochar='%s_%s'%(iz+1,jz+1)
        ## for 309c
        Csysfile='%s_%s_%s.dat'%(filetop_309c,tomochar,LFVER_309c)
        Csysdata = np.loadtxt(Csysfile)
        # used columns
        theta=Csysdata[:,1]
        Csys_p=Csysdata[:,3]
        err_Csys_p=Csysdata[:,5]
        epsfepsf=Csysdata[:,9]
        gepsf=Csysdata[:,12]
        npairs_weighted=Csysdata[:,20]
        del Csysdata
        #Lets see if we can't do something better for the CSys errors, as the Treecorr estimate doesn't
        #take the weights into account
        sigma_epsf=0.021667148635665007
        sigma_e=0.384   #both components
        err_gepsf=sigma_e*sigma_epsf/np.sqrt(npairs_weighted)
        err_epsfepsf=sigma_epsf*sigma_epsf/np.sqrt(npairs_weighted)
        an_err_Csys_p = Csys_p * np.sqrt(2*(err_gepsf/gepsf)**2 + (err_epsfepsf/epsfepsf)**2)
        # and plot the results with annotations of the bin combination and p-value
        tmp = ax.errorbar(theta, Csys_p/xiptheory, yerr=an_err_Csys_p/xiptheory, 
            color='gray', linestyle='--')
        if gridpos == 0:
            handles.append(tmp)
            LABELs.append('KiDS-1000-v1')
        ## for 321
        Csysfile='%s_%s_%s.dat'%(filetop_321,tomochar,LFVER_321)
        Csysdata = np.loadtxt(Csysfile)
        # used columns
        theta=Csysdata[:,1]
        Csys_p=Csysdata[:,3]
        err_Csys_p=Csysdata[:,5]
        epsfepsf=Csysdata[:,9]
        gepsf=Csysdata[:,12]
        npairs_weighted=Csysdata[:,20]
        del Csysdata
        #Lets see if we can't do something better for the CSys errors, as the Treecorr estimate doesn't
        #take the weights into account
        sigma_epsf=0.021667148635665007
        sigma_e=0.384   #both components
        err_gepsf=sigma_e*sigma_epsf/np.sqrt(npairs_weighted)
        err_epsfepsf=sigma_epsf*sigma_epsf/np.sqrt(npairs_weighted)
        an_err_Csys_p = Csys_p * np.sqrt(2*(err_gepsf/gepsf)**2 + (err_epsfepsf/epsfepsf)**2)
        # and plot the results with annotations of the bin combination and p-value
        tmp = ax.errorbar(theta, Csys_p/xiptheory, yerr=an_err_Csys_p/xiptheory, 
            color='red')
        if gridpos == 0:
            handles.append(tmp)
            LABELs.append('KiDS-1000-v2')

        ax.annotate(tomochar, xy=(0.19,0.1),xycoords='axes fraction',
                    size=14, ha='right', va='top')
        ax.set_xscale('log')        
        if iz<1 or jz<2:
            ax.set_ylim(-0.3,0.35)
        elif iz==1:
            ax.set_ylim(-0.12,0.15)
        else:
            ax.set_ylim(-0.09,0.09)
        ax.set_xlim(0.5,300.0)

        ax.set_xticks([1, 10, 100]) 
        ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
        # only label the subplots at the edges
        ax.label_outer()

#add labels
fig.legend(handles, LABELs,
        loc = 'center', ncol=3,
        bbox_to_anchor=(0.55, 0.98), fancybox=True, shadow=True)

axes[2,0].set_ylabel(r'$\xi_+^{\rm sys}/\xi_+^{\Lambda{\rm CDM}}$',fontsize=21)

axes[4,0].set_xlabel(r'$\theta$ (arcmin)',fontsize=17)
axes[4,1].set_xlabel(r'$\theta$ (arcmin)',fontsize=17)
axes[4,2].set_xlabel(r'$\theta$ (arcmin)',fontsize=17)

plt.tight_layout()
fig.subplots_adjust(top=0.95)

if outpath == 'show':
    plt.show()
    plt.close()
else:
    plt.savefig(outpath, dpi=300)
    plt.close()
    print("plot saved as", outpath)