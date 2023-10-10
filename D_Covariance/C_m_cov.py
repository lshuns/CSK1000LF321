# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-04-18 18:45:43
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-24 11:59:45

### calculate the m bias uncertainties 
###### only consider the statistical uncertainties

import numpy as np
import pandas as pd

# the shear bias
inpath_m = '../B_Shear_Bias/results/m_skills_v07D7p1_LF_321_kidsPhotometry_A12_gold_Rewei_varCorr_dz0p1_2D_PSFmodellingCorr.csv'
cata_tmp = pd.read_csv(inpath_m)
## the first is for whole
sigma_m = (cata_tmp.loc[1:, 'm1_err'].values + cata_tmp.loc[1:, 'm2_err'].values)/2.
del cata_tmp

# where to save
outpath = '/disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/COV/input/m_cov_statistical_only.ascii'

# construct m cov matrix for input m uncertainties
m_cov = np.diag(sigma_m**2)

# output matrix
np.savetxt(outpath, m_cov, fmt='%.4e')
print('saved to', outpath)
