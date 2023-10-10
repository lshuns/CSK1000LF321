#!/usr/bin/bash

# @Author: lshuns
# @Date:   2023-04-27 21:21:45
# @Last Modified by:   lshuns
# @Last Modified time: 2023-06-26 10:38:21

# activate kcap
source /disks/shear10/ssli/K1000CS/my_kcap/env/bin/activate
source cosmosis-configure

# # fiducial with scale cuts
# cosmosis ./config/COSEBIs_covariance_n5_t2.ini

# # without scale cuts
# cosmosis ./config/COSEBIs_covariance_n5_t05.ini

# # for Csys check
# cosmosis ./config/COSEBIs_covariance_n20_t05.ini

# for Csys check with scale cut
cosmosis ./config/COSEBIs_covariance_n20_t2.ini