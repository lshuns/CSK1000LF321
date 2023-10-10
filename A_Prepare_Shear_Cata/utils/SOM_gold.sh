#!/usr/bin/env bash
#
# Script for running SOM gold selection
#   hacked from CosmoWrapper: https://github.com/lshuns/CosmoWrapper

set -e

# Source the default parameters file
SOURCEFILE=$1
echo "use parameter file ${SOURCEFILE}"
source ${SOURCEFILE}
export OMP_NUM_THREADS=${MAXTHREADS}

################################### STEP 01 ###################################
# Prepare the intput data: replace the original magnitudes by the adapted ones,
# convert to LDAC table format if necessary
photcat="${PHOTCAT_ALL%.*}.cat"
SPECCAT_ALL_ADAPT=${SPECCAT_ALL//.csv/_adapt.fits}
# Construct the Spectrosopic Adapt Catalogue
if [ ! -f ${ROOTDIR}/${OUTPUTDIR}/${SPECCAT_ALL_ADAPT} ]
then
    echo "Constructing Spectroscopic Adapt Catalogue {"
    ${P_RSCRIPT} ${ROOTDIR}/construct_adapt_catalogue.R ${ROOTDIR}/${INPUTDIR}/${SPECCAT_ALL} ${ROOTDIR}/${OUTPUTDIR}/${SPECCAT_ALL_ADAPT}
    echo "} - Done"
else
    echo "Spectroscopic Adapt Catalogue Already Exists! Skipping!"
fi
# Convert to LDAC if FITS
if [ ! -f ${ROOTDIR}/${OUTPUTDIR}/${photcat} ]
then
    has_fields=$(echo $(${DIR_LDAC}/ldacdesc${THELI} -i ${ROOTDIR}/${INPUTDIR}/${PHOTCAT_ALL} | grep -c FIELDS))
    if [ "${has_fields}" == "1" ]
    then
        # Link file
        echo "Link the LDAC photcat {"
        ln -sv ${ROOTDIR}/${INPUTDIR}/${PHOTCAT_ALL} ${ROOTDIR}/${OUTPUTDIR}/${photcat}
        echo "} - Done"
    else
        # Create a fake LDAC table
        echo "Convert FITS file to LDAC {"
        python ${ROOTDIR}/fits2ldac.py ${ROOTDIR}/${INPUTDIR}/${PHOTCAT_ALL} -o ${ROOTDIR}/${OUTPUTDIR}/${photcat}
        echo "} - Done"
    fi
else
    echo "Photometric Catalogue Already Exists! Skipping!"
fi

# substitute to expected file extension .cat
PHOTCAT_ALL=$photcat

################################### STEP 02 ###################################
# Keep a minimum set of columns in the shear catalogue
PHOTCAT_ALL_DCOL=${PHOTCAT_ALL//.cat/_DIRcols.cat}
# Select DIRcol subset (within LDAC)
if [ ! -f ${ROOTDIR}/${OUTPUTDIR}/${PHOTCAT_ALL_DCOL} ]
then
    echo -e "Constructing the DIR Column Photometry Catalogue {"
    list=`${DIR_LDAC}/ldacdesc${THELI} -i ${ROOTDIR}/${OUTPUTDIR}/${PHOTCAT_ALL} -t OBJECTS |
        grep "Key name" | awk -F. '{print $NF}' |
        grep -v "SeqNr\|THELI_\|_B\|ID\|MAG_GAAP_\|${WEIGHTNAME}\|MAG_AUTO"`
    ${DIR_LDAC}/ldacdelkey${THELI} -i ${ROOTDIR}/${OUTPUTDIR}/${PHOTCAT_ALL} -k ${list} -o ${ROOTDIR}/${OUTPUTDIR}/${PHOTCAT_ALL_DCOL}
    echo "} - Done"
else
    echo "DIR Column Photometry Catalogue Already Exists! Skipping!"
fi

################################### STEP 03 ###################################
# Train the SOM on the calibration data
## Construct the Fiducial SOM
if [ ! -f ${ROOTDIR}/${OUTPUTDIR}/${SOMFILE} ]
then
    echo "Constructing the Fiducial SOM {"
    ${P_RSCRIPT} ${ROOTDIR}/SOM_DIR.R \
        -r ${ROOTDIR}/${OUTPUTDIR}/${PHOTCAT_ALL_DCOL} -t ${ROOTDIR}/${OUTPUTDIR}/${SPECCAT_ALL_ADAPT} \
        --toroidal --topo hexagonal --som.dim 101 101 -np -fn Inf \
        -sc ${MAXTHREADS} --only.som \
        -o ${ROOTDIR}/${OUTPUTDIR} -of ${SOMFILE//_SOMdata/} \
        --zr.label Z_B --zt.label ${ZLABEL} \
        -k MAG_GAAP_u-MAG_GAAP_g \
        MAG_GAAP_u-MAG_GAAP_r MAG_GAAP_g-MAG_GAAP_r \
        MAG_GAAP_u-MAG_GAAP_i MAG_GAAP_g-MAG_GAAP_i \
        MAG_GAAP_r-MAG_GAAP_i MAG_GAAP_u-MAG_GAAP_Z \
        MAG_GAAP_g-MAG_GAAP_Z MAG_GAAP_r-MAG_GAAP_Z \
        MAG_GAAP_i-MAG_GAAP_Z MAG_GAAP_u-MAG_GAAP_Y \
        MAG_GAAP_g-MAG_GAAP_Y MAG_GAAP_r-MAG_GAAP_Y \
        MAG_GAAP_i-MAG_GAAP_Y MAG_GAAP_Z-MAG_GAAP_Y \
        MAG_GAAP_u-MAG_GAAP_J MAG_GAAP_g-MAG_GAAP_J \
        MAG_GAAP_r-MAG_GAAP_J MAG_GAAP_i-MAG_GAAP_J \
        MAG_GAAP_Z-MAG_GAAP_J MAG_GAAP_Y-MAG_GAAP_J \
        MAG_GAAP_u-MAG_GAAP_H MAG_GAAP_g-MAG_GAAP_H \
        MAG_GAAP_r-MAG_GAAP_H MAG_GAAP_i-MAG_GAAP_H \
        MAG_GAAP_Z-MAG_GAAP_H MAG_GAAP_Y-MAG_GAAP_H \
        MAG_GAAP_J-MAG_GAAP_H MAG_GAAP_u-MAG_GAAP_Ks \
        MAG_GAAP_g-MAG_GAAP_Ks MAG_GAAP_r-MAG_GAAP_Ks \
        MAG_GAAP_i-MAG_GAAP_Ks MAG_GAAP_Z-MAG_GAAP_Ks \
        MAG_GAAP_Y-MAG_GAAP_Ks MAG_GAAP_J-MAG_GAAP_Ks \
        MAG_GAAP_H-MAG_GAAP_Ks MAG_AUTO
    echo "} - Done"
else
    echo "Fiducial SOM Already Exists! Skipping!"
fi

################################### STEP 04 ###################################
# Define the gold sample and compute the SOM redshift distributions
# Construct the Gold Classes
if [ ! -f ${ROOTDIR}/${OUTPUTDIR}/${PHOTCAT_ALL_DCOL//.cat/_allgoldclass.fits} ]
then
    echo "Constructing the Goldclass subsets {"
    time ${P_RSCRIPT} ${ROOTDIR}/construct_dr4_goldclasses.R \
        -p ${ROOTDIR}/${OUTPUTDIR}/${PHOTCAT_ALL_DCOL} \
        -s ${ROOTDIR}/${OUTPUTDIR}/${SPECCAT_ALL_ADAPT} \
        --som ${ROOTDIR}/${OUTPUTDIR}/${SOMFILE} \
        --blinds ${BLINDS} \
        --outputpath ${ROOTDIR}/${OUTPUTDIR}/ \
        --nzformat .asc \
        -cv ${WEIGHTNAME} \
        -zn ${ZLABEL} 
    echo "} - Done"
else
    echo "Gold Catalogues Already Exists! Skipping!"
    echo "(${ROOTDIR}/${OUTPUTDIR}/${PHOTCAT_ALL_DCOL//.cat/_allgoldclass.fits})"
fi

################################### STEP 05 ###################################
# Merge all SOM weights and gold flags into one catalogue
PHOTCAT_ALL_GOLD=${PHOTCAT_ALL//.cat/_goldclasses.cat}
# Merge the new GoldClasses back with the original catalogue
if [ ! -f ${ROOTDIR}/${OUTPUTDIR}/${PHOTCAT_ALL_GOLD} ]
then
    echo "Constructing merged catalogue {"
    GOLDFLAGLIST=""
    for goldset in ${GOLDLIST}
    do
        GOLDFLAGLIST=`echo $GOLDFLAGLIST Flag_SOM_${goldset}`
    done
    python ${ROOTDIR}/merge_ldac_and_fits.py ${ROOTDIR}/${OUTPUTDIR}/${PHOTCAT_ALL} \
        ${ROOTDIR}/${OUTPUTDIR}/${PHOTCAT_ALL_DCOL//.cat/_allgoldclass.fits} \
        ${ROOTDIR}/${OUTPUTDIR}/${PHOTCAT_ALL_GOLD} "${BLINDS}" "${GOLDFLAGLIST}"
    echo "} - Done"
else
    echo "Merged catalogue already exists! Skipping!"
    echo "(${ROOTDIR}/${OUTPUTDIR}/${PHOTCAT_ALL_GOLD})"
fi