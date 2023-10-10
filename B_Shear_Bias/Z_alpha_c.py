# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-02-28 18:15:34
# @Last Modified by:   lshuns
# @Last Modified time: 2023-02-28 18:35:55

### calculate alpha and c (PSF leakage)

import numpy as np 
import pandas as pd 

from astropy.io import fits
import statsmodels.api as sm 

# >>>>>>>>>> the catalogue

# ## ++++++++ public 309c cata
# inpath = '/disks/shear15/KiDS/KiDS-1000/K1000-SHEAR-9bandPZ-CATALOGUES/KiDS_DR4.1_ugriZYJHKs_SOM_gold_WL_cat.fits'
# with fits.open(inpath) as hdul:
#     cata = hdul['OBJECTS'].data
# ### used columns
# cata = pd.DataFrame({'e1_out': cata['e1'].astype(float),
#                         'e2_out': cata['e2'].astype(float),
#                         'shape_weight': cata['weight'].astype(float),
#                         'e1_psf': cata['PSF_e1'].astype(float),
#                         'e2_psf': cata['PSF_e2'].astype(float),
#                         'Z_B': cata['Z_B'].astype(float)})
# ### where to save
# outpath = './results/alpha_K1000_LF309c_public.csv'
# # Elapsed:3:41.69,User=97.689,System=116.632,CPU=96.6%.
# ## ++++++++ public 309c cata

## ++++++++ the 321 cata
inpath = '/disks/shear10/ssli/K1000CS/LF321_Outputs/K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A1_goldclasses.cat.A2.feather'
cata = pd.read_feather(inpath)
### used columns
cata = pd.DataFrame({'e1_out': cata['AlphaRecalD2_e1'].astype(float),
                        'e2_out': cata['AlphaRecalD2_e2'].astype(float),
                        'shape_weight': cata['AlphaRecalC_weight'].astype(float),
                        'e1_psf': cata['PSF_e1'].astype(float),
                        'e2_psf': cata['PSF_e2'].astype(float),
                        'Z_B': cata['Z_B'].astype(float)})
### where to save
outpath = './results/alpha_K1000_LF321_alphaRecal.csv'
# Elapsed:1:13.05,User=131.598,System=188.077,CPU=437.5%.
## ++++++++ the 321 cata

# >>>>>>>>>> genearl setup
ZBbins = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2]

# >>>>>>>>>> the estimate function
def alphaCalFunc_least_squares(cataSim):
    """
    Calculate the alpha term for a given simulated catalogue
        by requiring sum(e_out - alpha*e_psf) -> min

        Used columns and their names:
            e1_out, e2_out: measured ellipticity
            e1_psf, e2_psf: measured PSF ellipticity
            shape_weight: shape measurement weights
    """

    # out shear
    e1_out = np.array(cataSim['e1_out'])
    e2_out = np.array(cataSim['e2_out'])
    wgSim = np.array(cataSim['shape_weight'])

    # out PSF 
    e1_psf = np.array(cataSim['e1_psf'])
    e2_psf = np.array(cataSim['e2_psf'])

    # get least square values
    ## e1
    mod_wls = sm.WLS(e1_out, sm.add_constant(e1_psf), weights=wgSim)
    res_wls = mod_wls.fit()
    alpha1 = res_wls.params[1]
    c1 = res_wls.params[0]
    alpha1_err = (res_wls.cov_params()[1, 1])**0.5
    c1_err = (res_wls.cov_params()[0, 0])**0.5
    ## e2
    mod_wls = sm.WLS(e2_out, sm.add_constant(e2_psf), weights=wgSim)
    res_wls = mod_wls.fit()
    alpha2 = res_wls.params[1]
    c2 = res_wls.params[0]
    alpha2_err = (res_wls.cov_params()[1, 1])**0.5
    c2_err = (res_wls.cov_params()[0, 0])**0.5

    # save
    res = {'alpha1': alpha1, 'alpha2': alpha2,
            'c1': c1, 'c2': c2,
            'alpha1_err': alpha1_err, 'alpha2_err': alpha2_err,
            'c1_err': c1_err, 'c2_err': c2_err,
            }

    return res

# >>>>>>>>> workhorse

# output
f = open(outpath, 'w')

# the whole results
## mean and median of the x axis
mean_bin = np.average(cata['Z_B'].values, weights=cata['shape_weight'].values)
## the whole results
res = alphaCalFunc_least_squares(cata)
N_s = len(cata)

# collect columns names and values
cols = ','.join(list(res.keys()))
vals = ','.join(["{0:0.4f}".format(val) for val in res.values()])
cols = cols + f',Z_B_min,Z_B_max,Z_B_mean,Nobj'
vals = vals + f',-999,999,{mean_bin},{N_s}'    
print(cols, file=f)
print(vals, file=f)

for i_bin in range(len(ZBbins)-1):
    min_bin = ZBbins[i_bin]
    max_bin = ZBbins[i_bin+1]
    mask_bin = (cata['Z_B']>min_bin) & (cata['Z_B']<=max_bin)

    cata_selec = cata[mask_bin].copy()
    del mask_bin
    cata_selec.reset_index(drop=True, inplace=True)

    # mean and median of the x axis
    mean_bin = np.average(cata_selec['Z_B'].values, weights=cata_selec['shape_weight'].values)

    res = alphaCalFunc_least_squares(cata_selec)
    N_s = len(cata_selec)
    del cata_selec

    # print out values
    vals = ','.join(["{0:0.4f}".format(val) for val in res.values()]) \
    + f',{min_bin},{max_bin},{mean_bin},{N_s}'
    print(vals, file=f)

f.close()
print('results saved to', outpath)