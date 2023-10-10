#!/usr/bin/bash

# @Author: lshuns
# @Date:   2023-01-05 17:17:15
# @Last Modified by:   lshuns
# @Last Modified time: 2023-02-27 10:12:51

#### calculate the raw m (without correction of the 'shear interplay' effect and PSF modelling errors) 
####        with gold selection
####        in tomographic bins
####        without data reweighting
#### require biasEstimation code, available from https://github.com/KiDS-WL/MultiBand_ImSim/tree/main/biasEstimation

python /net/grecht/data2/ssli_files/Projects/Projects/8ShearBias_ImSim/MultiBand_ImSim/biasEstimation/bias_estimate_script.py\
        --bias_type m_LS_tile\
        --in_file /disks/shear10/ssli/K1000CS/skills_v07D7p1_Outputs/skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_newCut_A1_WeiCut_goldclasses.cat.A2.feather\
        --out_path ./results/m_skills_v07D7p1_LF_321_kidsPhotometry_A12_gold.csv\
        --cols_e12 AlphaRecalD2_e1 AlphaRecalD2_e2\
        --e_type measured\
        --cols_g12 g1_in g2_in\
        --col_weight AlphaRecalC_weight\
        --col_label tile_label\
        --col_binning Z_B\
        --binning_edges 0.1 0.3 0.5 0.7 0.9 1.2