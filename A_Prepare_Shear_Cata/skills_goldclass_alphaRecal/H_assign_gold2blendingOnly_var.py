# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-01-20 14:40:16
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-22 17:21:41

### assign the gold flag to varShear catalogues

import numpy as np
import pandas as pd 

# >>>>>>>>>>>> general info

# the constShear catalogue with gold selection
inpath_const = '/disks/shear10/ssli/K1000CS//skills_v07D7p1_Outputs//skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_newCut_A1_WeiCut_goldclasses.cat.A2.feather'
# the varShear catalogues
inpath_var = '/disks/shear16/ssli/ImSim/output/skills_v07D7v/skills_v07D7p1v_LF_321_shear_noSG_noWeiCut_newCut_var_withextra.feather'

outpath = inpath_var.replace('.feather', '_goldSelected.feather')

unique_shear_tags = ['CS0_rot0',  'CS0_rot90', 'CS0_rot180', 'CS0_rot270', 'CS0_rot45',  'CS0_rot135', 'CS0_rot225', 'CS0_rot315',
                    'CS1_rot0',  'CS1_rot90', 'CS1_rot180', 'CS1_rot270', 'CS1_rot45',  'CS1_rot135', 'CS1_rot225', 'CS1_rot315',
                    'CS2_rot0',  'CS2_rot90', 'CS2_rot180', 'CS2_rot270', 'CS2_rot45',  'CS2_rot135', 'CS2_rot225', 'CS2_rot315',
                    'CS3_rot0',  'CS3_rot90', 'CS3_rot180', 'CS3_rot270', 'CS3_rot45',  'CS3_rot135', 'CS3_rot225', 'CS3_rot315']
unique_shear_tags_ori = ['m283p283', 'p283p283', 'p283m283', 'm283m283', 'm283p283', 'p283p283', 'p283m283', 'm283m283',
                        'm283p283', 'p283p283', 'p283m283', 'm283m283', 'm283p283', 'p283p283', 'p283m283', 'm283m283',
                        'm283p283', 'p283p283', 'p283m283', 'm283m283', 'm283p283', 'p283p283', 'p283m283', 'm283m283',
                        'm283p283', 'p283p283', 'p283m283', 'm283m283', 'm283p283', 'p283p283', 'p283m283', 'm283m283']

# >>>>>>>>>>>> workhorse

# load catalogues
cata_const = pd.read_feather(inpath_const)[['run_tag', 'tile_label', 'gal_rot', 'id_input']]
print('number constShear', len(cata_const))
cata_var = pd.read_feather(inpath_var)
print('number varShear', len(cata_var))

# select those only in constShear gold cata
Ntot = 0
cata_final = []
for i_shear, shear_tag in enumerate(unique_shear_tags):

    print('for tag', shear_tag)

    # select constShear
    cata0_selec = cata_const.loc[cata_const['run_tag'] == unique_shear_tags_ori[i_shear], ['tile_label', 'gal_rot', 'id_input']]
    print('>>> number constShear', len(cata0_selec))

    # select varShear
    cata_selec = cata_var[cata_var['run_tag'] == shear_tag]
    Nbefore = len(cata_selec)
    Ntot += Nbefore
    print('>>> number varShear', Nbefore)

    # get gold 
    cata_selec = cata_selec.merge(cata0_selec, on=['tile_label', 'gal_rot', 'id_input'])
    del cata0_selec
    print('>>> number gold', len(cata_selec))

    # collect
    cata_final.append(cata_selec)
    del cata_selec

print('total number original', Ntot)
cata_final = pd.concat(cata_final, ignore_index=True)
print(f'Total number gold {len(cata_final)}')
cata_final.to_feather(outpath)
print('combined cata saved to', outpath)

# (py377) amsteldiep [167] > python H_assign_gold2blendingOnly_var.py 
# number constShear 28885361
# number varShear 128229095
# for tag CS0_rot0
# >>> number constShear 7222023
# >>> number varShear 4004960
# >>> number gold 2391399
# for tag CS0_rot90
# >>> number constShear 7222365
# >>> number varShear 4009125
# >>> number gold 2393830
# for tag CS0_rot180
# >>> number constShear 7221808
# >>> number varShear 4008918
# >>> number gold 2393857
# for tag CS0_rot270
# >>> number constShear 7219165
# >>> number varShear 4006144
# >>> number gold 2393912
# for tag CS0_rot45
# >>> number constShear 7222023
# >>> number varShear 4005099
# >>> number gold 2393031
# for tag CS0_rot135
# >>> number constShear 7222365
# >>> number varShear 4008387
# >>> number gold 2392268
# for tag CS0_rot225
# >>> number constShear 7221808
# >>> number varShear 4008826
# >>> number gold 2395374
# for tag CS0_rot315
# >>> number constShear 7219165
# >>> number varShear 4005494
# >>> number gold 2392094
# for tag CS1_rot0
# >>> number constShear 7222023
# >>> number varShear 4005800
# >>> number gold 2393864
# for tag CS1_rot90
# >>> number constShear 7222365
# >>> number varShear 4008517
# >>> number gold 2391780
# for tag CS1_rot180
# >>> number constShear 7221808
# >>> number varShear 4009210
# >>> number gold 2396320
# for tag CS1_rot270
# >>> number constShear 7219165
# >>> number varShear 4005964
# >>> number gold 2391751
# for tag CS1_rot45
# >>> number constShear 7222023
# >>> number varShear 4005113
# >>> number gold 2392970
# for tag CS1_rot135
# >>> number constShear 7222365
# >>> number varShear 4008391
# >>> number gold 2392182
# for tag CS1_rot225
# >>> number constShear 7221808
# >>> number varShear 4009000
# >>> number gold 2395269
# for tag CS1_rot315
# >>> number constShear 7219165
# >>> number varShear 4006345
# >>> number gold 2392217
# for tag CS2_rot0
# >>> number constShear 7222023
# >>> number varShear 4005013
# >>> number gold 2391431
# for tag CS2_rot90
# >>> number constShear 7222365
# >>> number varShear 4009464
# >>> number gold 2394158
# for tag CS2_rot180
# >>> number constShear 7221808
# >>> number varShear 4008662
# >>> number gold 2393631
# for tag CS2_rot270
# >>> number constShear 7219165
# >>> number varShear 4006548
# >>> number gold 2393906
# for tag CS2_rot45
# >>> number constShear 7222023
# >>> number varShear 4003785
# >>> number gold 2389546
# for tag CS2_rot135
# >>> number constShear 7222365
# >>> number varShear 4009269
# >>> number gold 2395790
# for tag CS2_rot225
# >>> number constShear 7221808
# >>> number varShear 4007872
# >>> number gold 2391694
# for tag CS2_rot315
# >>> number constShear 7219165
# >>> number varShear 4006773
# >>> number gold 2395593
# for tag CS3_rot0
# >>> number constShear 7222023
# >>> number varShear 4004071
# >>> number gold 2388998
# for tag CS3_rot90
# >>> number constShear 7222365
# >>> number varShear 4009608
# >>> number gold 2396577
# for tag CS3_rot180
# >>> number constShear 7221808
# >>> number varShear 4008624
# >>> number gold 2391529
# for tag CS3_rot270
# >>> number constShear 7219165
# >>> number varShear 4007017
# >>> number gold 2396465
# for tag CS3_rot45
# >>> number constShear 7222023
# >>> number varShear 4003823
# >>> number gold 2389569
# for tag CS3_rot135
# >>> number constShear 7222365
# >>> number varShear 4009206
# >>> number gold 2395692
# for tag CS3_rot225
# >>> number constShear 7221808
# >>> number varShear 4007679
# >>> number gold 2391669
# for tag CS3_rot315
# >>> number constShear 7219165
# >>> number varShear 4006388
# >>> number gold 2395581
# total number original 128229095
# Total number gold 76583947
# combined cata saved to /disks/shear16/ssli/ImSim/output/skills_v07D7v/skills_v07D7p1v_LF_321_shear_noSG_noWeiCut_newCut_var_withextra_goldSelected.feather
# Elapsed:23:58.21,User=1344.103,System=2480.875,CPU=265.9%.
