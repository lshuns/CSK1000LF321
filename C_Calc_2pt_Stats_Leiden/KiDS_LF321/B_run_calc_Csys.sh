# @Author: lshuns
# @Date:   2023-02-23 09:21:35
# @Last Modified by:   lshuns
# @Last Modified time: 2023-03-09 17:35:39

### run the calc Csys using Catherine's code
for ZBIN1 in `seq 5`
do      
    for ZBIN2 in `seq $ZBIN1 5`
    do
        echo "working on" $ZBIN1 $ZBIN2
        python calc_Csys_stargal_w_treecorr.py 9 0.5 300 \
                /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/TOMOCATS/K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_5Z_${ZBIN1}.fits\
                /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/TOMOCATS/K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_5Z_${ZBIN1}.fits\
                /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/TOMOCATS/K1000_N_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_5Z_${ZBIN2}.fits\
                /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/TOMOCATS/K1000_S_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_5Z_${ZBIN2}.fits\
                /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/Csys/CSys_BLIND_NONE_5Z_${ZBIN1}_${ZBIN2}_LF_glab_321_v2_A12_goldclasses.dat
    done
done
# Elapsed:33:58.48,User=36101.603,System=1773.593,CPU=1858.0%.
