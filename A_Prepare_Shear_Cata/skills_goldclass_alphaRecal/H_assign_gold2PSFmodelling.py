# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-04-11 17:31:44
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-22 17:25:09

### assign the gold flag to PSF modelling test catalogues

import numpy as np
import pandas as pd 

# >>>>>>>>>>>> general info

# the constShear catalogue with gold selection
inpath0 = '/disks/shear10/ssli/K1000CS//skills_v07D7p1_Outputs//skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_newCut_A1_WeiCut_goldclasses.cat.A2.feather'
# the PSF modelling test catalogues
inpath_test = '/disks/shear16/ssli/ImSim/output/skills_v07D7p1_PSFmodelling/skills_v07D7p1_LF_321_combined_PSFmodelling41_321_noSG_noWeiCut_newCut.feather'

outpath = inpath_test.replace('.feather', '_goldSelected.feather')

# >>>>>>>>>>>> workhorse

# load catalogues
cata0 = pd.read_feather(inpath0)[['run_tag', 'tile_label', 'gal_rot', 'NUMBER']]
print('number whole', len(cata0))
cata_test = pd.read_feather(inpath_test)
print('number PSFmodelling', len(cata_test))

# select those only in constShear gold cata
cata_test = cata_test.merge(cata0, on=['run_tag', 'tile_label', 'gal_rot', 'NUMBER'])
del cata0
print('>>> number gold', len(cata_test))

# save
cata_test.to_feather(outpath)
print('combined cata saved to', outpath)

# (py377) amsteldiep [168] > python H_assign_gold2PSFmodelling.py 
# number whole 28885361
# number PSFmodelling 13787117
# >>> number gold 8291051
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7p1_PSFmodelling/skills_v07D7p1_LF_321_combined_PSFmodelling41_321_noSG_noWeiCut_newCut_goldSelected.feather
# Elapsed:2:22.45,User=168.439,System=655.451,CPU=578.3%.
