# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-01-20 17:23:20
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-22 16:54:35

### assign the gold flag to blending-only constShear catalogues

import numpy as np
import pandas as pd 

# >>>>>>>>>>>> general info

# the constShear catalogue with gold selection
inpath0 = '/disks/shear10/ssli/K1000CS//skills_v07D7p1_Outputs//skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_newCut_A1_WeiCut_goldclasses.cat.A2.feather'
# the blending-only constShear catalogues
inpath_b = '/disks/shear16/ssli/ImSim/output/skills_v07D7v/skills_v07D7p1v_LF_321_shear_noSG_noWeiCut_newCut_const.feather'

outpath = inpath_b.replace('.feather', '_goldSelected.feather')

# >>>>>>>>>>>> workhorse

# load catalogues
cata0 = pd.read_feather(inpath0)[['run_tag', 'tile_label', 'gal_rot', 'id_input']]
print('number whole', len(cata0))
cata_b = pd.read_feather(inpath_b)
print('number blending-only', len(cata_b))

# select those only in constShear gold cata
cata_b = cata_b.merge(cata0, on=['run_tag', 'tile_label', 'gal_rot', 'id_input'])
del cata0
print('>>> number gold', len(cata_b))

# save
cata_b.to_feather(outpath)
print('combined cata saved to', outpath)

# (py377) amsteldiep [166] > python H_assign_gold2blendingOnly_const.py 
# number whole 28885361
# number blending-only 16025052
# >>> number gold 9593406
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7v/skills_v07D7p1v_LF_321_shear_noSG_noWeiCut_newCut_const_goldSelected.feather
# Elapsed:3:34.54,User=161.223,System=587.432,CPU=348.9%.
