# @Author: lshuns
# @Date:   2023-04-15 14:01:12
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-15 12:21:14

### combine the two-point correlation functions
###### This is a hack from https://github.com/KiDS-WL/Cat_to_Obs_K1000_P1/blob/master/Calc_2pt_Stats/doall_calc2pt.sh
######### mode: \"COMBINEXI\": combine the xi_+/- results from N/S for tomo bin pair i j"

# >>>>>>>>>>>> I/O

# input info
INDIR="/disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI"
# FILETAG="BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins"
# FILETAG="BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_2.0_300.0_zbins"
FILETAG="BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_9_theta_0.5_300.0_zbins"

# number of ZB bins
NZBINS=5

# >>>>>>>>>>>>> workhorse

# make a tmpdir for intermediate outputs
TMPDIR="${INDIR}/tmp"
mkdir ${TMPDIR}

for IZBIN in `seq ${NZBINS}`
do      
    for JZBIN in `seq ${IZBIN} ${NZBINS}`
    do

    echo "for tomo bin pair ${IZBIN} ${JZBIN}"

    ## Define paths
    inFileN="${INDIR}/XI_K1000_N_${FILETAG}_${IZBIN}_${JZBIN}.asc"
    inFileS="${INDIR}/XI_K1000_S_${FILETAG}_${IZBIN}_${JZBIN}.asc"
    outPath="${INDIR}/XI_K1000_ALL_${FILETAG}_${IZBIN}_${JZBIN}.asc"

    # check do the files exist?
    test -f ${inFileN} || \
    { echo "Error: KiDS-N XI results ${inFileN} do not exist!"; exit 1; } 
    test -f ${inFileS} || \
    { echo "Error: KiDS-S XI results ${inFileS} do not exist!"; exit 1; } 

    # and now lets combine them using the fabulous awk
    # which Tilman will be most scathing about, but I love it nevertheless
    # first lets grab the header which we want to replicate

    head -2 < ${inFileN} > ${TMPDIR}/xi_header

    # paste the two catalogues together
    paste ${inFileN} ${inFileS} > ${TMPDIR}/xi_paste

    # time for awk where we use npairs to weight every other
    # column to get the average
    # $10 = npairs in KiDS-N,  $20 = npairs in KiDS-S
    # For the sigma_xi column I'm assuming ngals in N and S are similar and sum sigma_xi in quadrature
    # This isn't correct but we don't really use the sigma_xi column ($8, $9, $19, and $20) 
    # Finally give the sum of weights and the sum of npairs
    
    awk 'NR>2 {printf " % .12e  % .12e  % .12e  % .12e  % .12e  % .12e  % .12e  % .12e  % .12e  % .12e  % .12e\n", 
                     ($1*$11 + $12*$22)/($11+$22), \
                     ($2*$11 + $13*$22)/($11+$22), \
                     ($3*$11 + $14*$22)/($11+$22), \
                     ($4*$11 + $15*$22)/($11+$22), \
                     ($5*$11 + $16*$22)/($11+$22), \
                     ($6*$11 + $17*$22)/($11+$22), \
                     ($7*$11 + $18*$22)/($11+$22), \
                     sqrt($8*$8 + $19*$19), \
                     sqrt($9*$9 + $20*$20), \
                     $10+$21, \
                     $11+$22}' < ${TMPDIR}/xi_paste > ${TMPDIR}/xi_comb
    
    # Finally put the header back
    cat ${TMPDIR}/xi_header ${TMPDIR}/xi_comb > ${outPath}

    # Did it work?
    test -f ${outPath} || \
        { echo "Error: Combined Treecorr output ${outPath} was not created! !"; exit 1; }

    echo "combined catalogue saved as ${outPath}"

    done
done

# remove tmp 
rm -rf ${TMPDIR}
echo "${TMPDIR} removed"

# for tomo bin pair 1 1
# combined catalogue saved as /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_ALL_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_1_1.asc
# for tomo bin pair 1 2
# combined catalogue saved as /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_ALL_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_1_2.asc
# for tomo bin pair 1 3
# combined catalogue saved as /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_ALL_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_1_3.asc
# for tomo bin pair 1 4
# combined catalogue saved as /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_ALL_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_1_4.asc
# for tomo bin pair 1 5
# combined catalogue saved as /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_ALL_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_1_5.asc
# for tomo bin pair 2 2
# combined catalogue saved as /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_ALL_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_2_2.asc
# for tomo bin pair 2 3
# combined catalogue saved as /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_ALL_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_2_3.asc
# for tomo bin pair 2 4
# combined catalogue saved as /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_ALL_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_2_4.asc
# for tomo bin pair 2 5
# combined catalogue saved as /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_ALL_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_2_5.asc
# for tomo bin pair 3 3
# combined catalogue saved as /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_ALL_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_3_3.asc
# for tomo bin pair 3 4
# combined catalogue saved as /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_ALL_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_3_4.asc
# for tomo bin pair 3 5
# combined catalogue saved as /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_ALL_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_3_5.asc
# for tomo bin pair 4 4
# combined catalogue saved as /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_ALL_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_4_4.asc
# for tomo bin pair 4 5
# combined catalogue saved as /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_ALL_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_4_5.asc
# for tomo bin pair 5 5
# combined catalogue saved as /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/XI_K1000_ALL_BLIND_NONE_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A12_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_5_5.asc
# /disks/shear10/ssli/K1000CS/LF321_Cat_to_Obs_Outputs/OUTSTATS/XI/tmp removed
