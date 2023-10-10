# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-05-06 13:36:47
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-21 17:04:03

### create necessary files for local minimisation
###### 1. pipeline.ini for each variation
###### 2. start point from fiducial minimisation run

# >>>>>>>>>>>> I/O

## the run name
# RUN_NAME = 'multinest_nIA'
RUN_NAME = 'polychord_nIA'
## variations for data vector
theta_tag = '2.00_300.00'
m_tag_list = ['no_m_bias', 'with_raw_m_bias',
                'test_m_nU', 'test_m_nD',
                'test_m_qU', 'test_m_qD',
                'test_m_sizeU', 'test_m_sizeD'
                ]
inpath_data_template ='/disks/shear10/ssli/K1000CS/LF321_cosmo/cosebis/data_vector/cosebis_K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_goldclasses_Flag_SOM_Fid_filt_{:}_nbins5_theta_{:}_cosmosis.fits'

## the original .ini file
inpath0 = f'../LF321_{RUN_NAME}/pipeline.ini'

## where to save
outpath_prefix = f'./config/{RUN_NAME}_pipeline'

# >>>>>>>>>>>>> workhorse

# read the file
with open(inpath0, 'r') as file:
    lines0 = file.readlines() 
# loop over variations
for m_tag in m_tag_list:
    inpath_data = inpath_data_template.format(m_tag, theta_tag)

    # modify the data_file
    lines = lines0.copy()
    for index, line in enumerate(lines):
        if 'data_file' in line:
            lines[index] = 'data_file = ' + inpath_data
            break

    # save the new ini file
    outpath = f'{outpath_prefix}_{m_tag}.ini'
    with open(outpath, 'w') as file:
        file.writelines(lines)
    print('new ini saved to ', outpath)