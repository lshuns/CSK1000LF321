# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-04-03 11:26:25
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-22 17:15:57

### calculate the shape noise and effective number density
#### reference: Joachimi et al. (2021): KiDS-1000 methodology
###### shape noise estimates (C.9)
###### effective number density (C.12)

import sys
import argparse

import numpy as np
import pandas as pd

from astropy.io import fits


# >>> general info

# ### 309c public
# inpath = '/disks/shear10/ssli/KiDS/K1000_blind_gold/K1000_combined_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_THELI_INT.fits'
# wei_col = 'recal_weight_C'
# e1_col = 'autocal_e1_C'
# e2_col = 'autocal_e2_C'
# gold_col = 'Flag_SOM_Fid_C'
# # number original 31446584
# # number gold 21262011
# # bin (0.1, 0.3]: 0.6163804051584554 0.2696215845705293
# # bin (0.3, 0.5]: 1.1821302809254923 0.2578848246158315
# # bin (0.5, 0.7]: 1.8541169742317374 0.27258637747974623
# # bin (0.7, 0.9]: 1.2592061444713252 0.25393698433475576
# # bin (0.9, 1.2]: 1.3108856755180338 0.27027412027375947
# # Elapsed:3:35.48,User=179.005,System=38.690,CPU=101.0%.

# ### 309c Jan Luca
# inpath = '/disks/shear10/ssli/K1000CS/LF309c_Outputs/K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses.cat'
# wei_col = 'recal_weight_C'
# e1_col = 'autocal_e1_C'
# e2_col = 'autocal_e2_C'
# gold_col = 'Flag_SOM_Fid_C'
# # number original 31446584
# # number gold 23427938
# # bin (0.1, 0.3]: 0.6873508919309291 0.2714744350433925
# # bin (0.3, 0.5]: 1.345926248029499 0.2620993559308912
# # bin (0.5, 0.7]: 1.997766285493571 0.27444731676935896
# # bin (0.7, 0.9]: 1.3745064357213232 0.2571866583427248
# # bin (0.9, 1.2]: 1.3208890533636062 0.2706390109489028
# # Elapsed:3:26.02,User=143.453,System=50.762,CPU=94.2%.

### 321
inpath = '/disks/shear10/ssli/K1000CS/LF321_Outputs//K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A1_goldclasses.cat.A2.feather'
wei_col = 'AlphaRecalC_weight'
e1_col = 'AlphaRecalD2_e1'
e2_col = 'AlphaRecalD2_e2'
gold_col = 'Flag_SOM_Fid_NONE'
# number original 23401764
# number gold 23401764
# bin (0.1, 0.3]: 0.6752667818128969 0.2702116434341894
# bin (0.3, 0.5]: 1.2979600496986519 0.26157689022690395
# bin (0.5, 0.7]: 1.9678242285798497 0.27651381922776336
# bin (0.7, 0.9]: 1.3931149418679294 0.2654044829986481
# bin (0.9, 1.2]: 1.347039366276972 0.28608453246947557
# Elapsed:1:06.80,User=82.673,System=706.549,CPU=1181.4%.

# effective area
Aeff = 777.4 * 60**2 # square arcmin

# binning col
bin_col = 'Z_B'
binning_edges = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2]

# >>> the function
def ShapeNoiseFunc(e1, e2, weight, 
                    sky_area=1, m_bias = 0,
                    c12=None):
    """
    calculate the effective number density
    """

    # the effective number density
    Neff = (np.sum(weight * (1+m_bias))**2 
            / np.sum(np.square(weight * (1+m_bias))) 
            / sky_area)
    # print('effective number density', Neff)

    # c bias
    if c12 is None:
        c1 = np.average(e1, weights=weight)
        c2 = np.average(e2, weights=weight)
        # print('calculated c1, c2', c1, c2)
    else:
        c1, c2 = c12
        # print('provided c1, c2', c1, c2)

    # e sigma
    wei2_withM = np.square(weight * (1+m_bias))
    sigma_e = (  # geometric mean
        np.sqrt(  # standard deviation of e1
          np.sum(np.square(weight * (e1 - c1))) /
          np.sum(wei2_withM)
        ) / 2.0 +
        np.sqrt(  # standard deviation of e2
          np.sum(np.square(weight * (e2 - c2))) /
          np.sum(wei2_withM)
        ) / 2.0)
    # print(f'sigma e', sigma_e)

    # return results
    res = {'Neff': Neff, 'sigma_e': sigma_e}

    return res

# >>> workhorse

## load cata
file_type = inpath[-3:]
if file_type == 'her':
    cata = pd.read_feather(inpath)
elif file_type == 'cat':
    with fits.open(inpath) as hdul:
        cata = hdul['OBJECTS'].data
elif file_type == 'its':
    with fits.open(inpath) as hdul:
        cata = hdul[1].data
else:
    raise Exception(f'Not supported input file type! {inpath}')
print('number original', len(cata))
### select gold
cata = cata[cata[gold_col]==1]
print('number gold', len(cata))
### used columns
cata = pd.DataFrame({'weight': np.array(cata[wei_col], dtype=float),
                        'e1': np.array(cata[e1_col], dtype=float),
                        'e2': np.array(cata[e2_col], dtype=float),
                        'bin_col': np.array(cata[bin_col], dtype=float)
                        })

## calculate for each bin
for i_bin in range(len(binning_edges)-1):
    bin_min = binning_edges[i_bin]
    bin_max = binning_edges[i_bin + 1]

    cata_selec = cata[(cata['bin_col']>bin_min)&(cata['bin_col']<=bin_max)]

    res = ShapeNoiseFunc(cata_selec['e1'].values, cata_selec['e2'].values, cata_selec['weight'].values, 
                    sky_area=Aeff, m_bias = 0,
                    c12=None)

    print(f'bin ({bin_min}, {bin_max}]:', res['Neff'], res['sigma_e'])
