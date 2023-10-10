#!/usr/bin/bash

# @Author: lshuns
# @Date:   2023-01-05 17:17:15
# @Last Modified by:   lshuns
# @Last Modified time: 2023-02-25 09:49:34

#### calculate the raw m (without correction of the 'shear interplay' effect and PSF modelling errors) 
####        with gold selection
####        in tomographic bins
####        reweighting the mock catalogue using data
#### require biasEstimation code, available from https://github.com/KiDS-WL/MultiBand_ImSim/tree/main/biasEstimation

python /net/grecht/data2/ssli_files/Projects/Projects/8ShearBias_ImSim/MultiBand_ImSim/biasEstimation/bias_estimate_script_dataRewei.py\
                    --bias_type m_LS_DataRewei_2D\
                    --in_file /disks/shear10/ssli/K1000CS/skills_v07D7p1_Outputs/skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_newCut_A1_WeiCut_goldclasses.cat.A2.feather\
                    --out_path ./results/m_skills_v07D7p1_LF_321_kidsPhotometry_A12_gold_Rewei.csv\
                    --cols_e12 AlphaRecalD2_e1 AlphaRecalD2_e2\
                    --cols_g12 g1_in g2_in\
                    --col_goldFlag_sim_data Flag_SOM_Fid_NONE Flag_SOM_Fid_NONE\
                    --col_weight_sim_data AlphaRecalC_weight AlphaRecalC_weight\
                    --col_label tile_label\
                    --col_binning_sim_data Z_B Z_B\
                    --binning_edges 0.1 0.3 0.5 0.7 0.9 1.2\
                    --in_file_data /disks/shear10/ssli/K1000CS/LF321_Outputs/K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A1_goldclasses.cat.A2.feather\
                    --bin1_info SNR_LF_r model_SNratio 20\
                    --bin2_info R R 20\
                    --fitting_method tile_based\
                    --save_surface_prefix ./results/m_skills_v07D7p1_LF_321_kidsPhotometry_A12_gold_Rewei\
                    --save_bounds_prefix ./results/m_skills_v07D7p1_LF_321_kidsPhotometry_A12_gold_Rewei