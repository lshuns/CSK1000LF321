# ----------------------------------------------------------------
# File Name:           calc_Csys_stargal_w_treecorr.py
# Author:              Catherine Heymans (heymans@roe.ac.uk)
# Description:         short python script to run treecorr to calculate CSysy
#                      given a KiDS fits catalogue
#                      script will need to change if keywords in KIDS cats are updated
#                      if I had more time, this would be far prettier.... :(
# ----------------------------------------------------------------

### Script history information:
# 24/02/2023 SSLi: new version of treecorr removed gen_write, copy that from old version

import treecorr
import sys
import numpy as np


# Read in user input to set the nbins, theta_min, theta_max, lin_not_log, fitscat1, fitscat2, outfilename
if len(sys.argv) <8: 
    print("Usage: %s nbins theta_min(arcmin) theta_max(arcmin) catalogueN_i.fits catalogueS_i.fits catalogueN_j.fits catalogueS_j.fits outfilename" % sys.argv[0]) 
    sys.exit(1)
else:
    nbins = int(sys.argv[1]) 
    theta_min = float(sys.argv[2]) 
    theta_max = float(sys.argv[3]) 
    fitscatN_i = sys.argv[4]
    fitscatS_i = sys.argv[5]
    fitscatN_j = sys.argv[6]
    fitscatS_j = sys.argv[7]
    outfile = sys.argv[8]

# prepare the catalogues

Ncat1_i = treecorr.Catalog(fitscatN_i, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  g1_col='e1', g2_col='e2', w_col='weight')
Ncat2_i = treecorr.Catalog(fitscatN_i, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  g1_col='PSF_e1', g2_col='PSF_e2', w_col='weight')
Ncat3_i = treecorr.Catalog(fitscatN_i, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  g1_col='dPSF_e1_xy', g2_col='dPSF_e2_xy', w_col='weight')
Ncat_sq_i = treecorr.Catalog(fitscatN_i, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  g1_col='e1', g2_col='e2', w_col='weightsq')

Scat1_i = treecorr.Catalog(fitscatS_i, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  g1_col='e1', g2_col='e2', w_col='weight')
Scat2_i = treecorr.Catalog(fitscatS_i, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  g1_col='PSF_e1', g2_col='PSF_e2', w_col='weight')
Scat3_i = treecorr.Catalog(fitscatS_i, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  g1_col='dPSF_e1_xy', g2_col='dPSF_e2_xy', w_col='weight')
Scat_sq_i = treecorr.Catalog(fitscatS_i, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  g1_col='e1', g2_col='e2', w_col='weightsq')

Ncat1_j = treecorr.Catalog(fitscatN_j, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  g1_col='e1', g2_col='e2', w_col='weight')
Ncat2_j = treecorr.Catalog(fitscatN_j, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  g1_col='PSF_e1', g2_col='PSF_e2', w_col='weight')
Ncat3_j = treecorr.Catalog(fitscatN_j, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  g1_col='dPSF_e1_xy', g2_col='dPSF_e2_xy', w_col='weight')
Ncat_sq_j = treecorr.Catalog(fitscatN_j, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  g1_col='e1', g2_col='e2', w_col='weightsq')

Scat1_j = treecorr.Catalog(fitscatS_j, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  g1_col='e1', g2_col='e2', w_col='weight')
Scat2_j = treecorr.Catalog(fitscatS_j, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  g1_col='PSF_e1', g2_col='PSF_e2', w_col='weight')
Scat3_j = treecorr.Catalog(fitscatS_j, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  g1_col='dPSF_e1_xy', g2_col='dPSF_e2_xy', w_col='weight')
Scat_sq_j = treecorr.Catalog(fitscatS_j, ra_col='ALPHA_J2000', dec_col='DELTA_J2000', ra_units='deg', dec_units='deg', \
                                  g1_col='e1', g2_col='e2', w_col='weightsq')

# when using fine bins this can increase to 1.5
inbinslop = 0.1

# Define the binning based on command line input
#star-gal #star-star # N and S
geN_ij = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', bin_slop=inbinslop)
geN_ji = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', bin_slop=inbinslop)
geS_ij = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', bin_slop=inbinslop)
geS_ji = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', bin_slop=inbinslop)

eeN = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', bin_slop=inbinslop)
eeS = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', bin_slop=inbinslop)

edeN = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', bin_slop=inbinslop)
edeS = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', bin_slop=inbinslop)

dedeN = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', bin_slop=inbinslop)
dedeS = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', bin_slop=inbinslop)

gdeN = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', bin_slop=inbinslop)
gdeS = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', bin_slop=inbinslop)

gewsqN = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', bin_slop=inbinslop)
gewsqS = treecorr.GGCorrelation(min_sep=theta_min, max_sep=theta_max, nbins=nbins, sep_units='arcmin', bin_slop=inbinslop)

# Calculate the 2pt correlation function
geN_ij.process(Ncat1_i,Ncat2_j)
geS_ij.process(Scat1_i,Scat2_j)

geN_ji.process(Ncat1_j,Ncat2_i)
geS_ji.process(Scat1_j,Scat2_i)

eeN.process(Ncat2_i,Ncat2_j)
eeS.process(Scat2_i,Scat2_j)

edeN.process(Ncat2_i,Ncat3_j)
edeS.process(Scat2_i,Scat3_j)

dedeN.process(Ncat3_i,Ncat3_j)
dedeS.process(Scat3_i,Scat3_j)

# just want this for the autobins at the moment
gdeN.process(Ncat1_i,Ncat3_i)
gdeS.process(Scat1_i,Scat3_i)


ge_p_all_ij = (geN_ij.xip*geN_ij.npairs + geS_ij.xip*geS_ij.npairs)/(geN_ij.npairs + geS_ij.npairs)
ge_m_all_ij = (geN_ij.xim*geN_ij.npairs + geS_ij.xim*geS_ij.npairs)/(geN_ij.npairs + geS_ij.npairs)
ge_p_all_ji = (geN_ji.xip*geN_ji.npairs + geS_ji.xip*geS_ji.npairs)/(geN_ji.npairs + geS_ji.npairs)
ge_m_all_ji = (geN_ji.xim*geN_ji.npairs + geS_ji.xim*geS_ji.npairs)/(geN_ji.npairs + geS_ji.npairs)

ee_p_all = (eeN.xip*eeN.npairs + eeS.xip*eeS.npairs)/(eeN.npairs + eeS.npairs)
ee_m_all = (eeN.xim*eeN.npairs + eeS.xim*eeS.npairs)/(eeN.npairs + eeS.npairs)

ede_p_all = (edeN.xip*edeN.npairs + edeS.xip*edeS.npairs)/(edeN.npairs + edeS.npairs)
ede_m_all = (edeN.xim*edeN.npairs + edeS.xim*edeS.npairs)/(edeN.npairs + edeS.npairs)

dede_p_all = (dedeN.xip*dedeN.npairs + dedeS.xip*dedeS.npairs)/(dedeN.npairs + dedeS.npairs)
dede_m_all = (dedeN.xim*dedeN.npairs + dedeS.xim*dedeS.npairs)/(dedeN.npairs + dedeS.npairs)

tot_npairs= gdeN.npairs + gdeS.npairs
gde_p_all = (gdeN.xip*gdeN.npairs + gdeS.xip*gdeS.npairs)/tot_npairs
gde_m_all = (gdeN.xim*gdeN.npairs + gdeS.xim*gdeS.npairs)/tot_npairs

sigma_gde_p_all = np.sqrt((gdeN.npairs/tot_npairs)**2 * gdeN.varxip + (gdeS.npairs/tot_npairs)**2 * gdeS.varxip)

#Combine to calculate Csys
Csys_p = (ge_p_all_ij*ge_p_all_ji)/(ee_p_all)
Csys_m = (ge_m_all_ij*ge_m_all_ji)/(ee_m_all)

#A shape noise estimate of the error
Csys_p_err = Csys_p * np.sqrt((geN_ij.varxip)/geN_ij.xip**2 + (geN_ji.varxip)/geN_ji.xip**2 + (eeN.varxip)/eeN.xip**2 + (geS_ij.varxip)/geS_ij.xip**2 + (geS_ji.varxip)/geS_ji.xip**2 + (eeS.varxip)/eeS.xip**2)
Csys_m_err = Csys_m * np.sqrt((geN_ij.varxim)/geN_ij.xim**2 + (geN_ij.varxim)/geN_ij.xip**2 + (eeN.varxim)/eeN.xim**2 + (geS_ij.varxim)/geS_ij.xim**2 + (geS_ji.varxim)/geS_ji.xim**2 + (eeS.varxim)/eeS.xim**2)

# we also need to calculate the weighted square 2pt correlation function to determine the effective number of pairs for analytical calculations
gewsqN.process(Ncat_sq_i,Ncat_sq_j)
gewsqS.process(Scat_sq_i,Scat_sq_j)

# Calculate the weighted Npairs = sum(weight_a*weight_b)^2 / sum(weight_a^2*weight_b^2)                                                                                                              
npairs_weighted = (geN_ij.weight)*(geN_ij.weight)/gewsqN.weight + (geS_ij.weight)*(geS_ij.weight)/gewsqS.weight 

#Use treecorr to write out the output file and praise-be once more for Jarvis and his well documented code
from treecorr.util import make_writer
def gen_write(file_name, col_names, columns, params=None, precision=4, file_type=None, logger=None):
    """Write some columns to an output file with the given column names.
    We do this basic functionality a lot, so put the code to do it in one place.
    :param file_name:   The name of the file to write to.
    :param col_names:   A list of columns names for the given columns.
    :param columns:     A list of numpy arrays with the data to write.
    :param params:      A dict of extra parameters to write at the top of the output file (for
                        ASCII output) or in the header (for FITS output).  (default: None)
    :param precision:   Output precision for ASCII. (default: 4)
    :param file_type:   Which kind of file to write to. (default: determine from the file_name
                        extension)
    :param logger:      If desired, a logger object for logging. (default: None)
    """
    if len(col_names) != len(columns):
        raise ValueError("col_names and columns are not the same length.")
    if len(columns) == 0:
        raise ValueError("len(columns) == 0")
    for col in columns[1:]:
        if col.shape != columns[0].shape:
            raise ValueError("columns are not all the same shape")

    columns = [ col.flatten() for col in columns ]
    writer = make_writer(file_name, precision, file_type, logger)
    with writer:
        writer.write(col_names, columns, params=params)

# treecorr.util.gen_write(outfile,
gen_write(outfile,
            ['r_nom','meanr','meanlogr','Csys_p','Csys_m', 'sigma_Csys_p','sigma_Csys_m','weight','npairs', 'xip_starstar', 'xim_starstar',
             'sigma_xi_starstar','xip_stargal_ij', 'xim_stargal_ij','xip_stargal_ji', 'xim_stargal_ji', 'xip_star_depsf', 'xip_depsf_depsf', 'xip_gal_depsf', 'sigma_xip_gal_depsf', 'npairs_weighted'],
            [ geN_ij.rnom,geN_ij.meanr, geN_ij.meanlogr,Csys_p,Csys_m,Csys_p_err,Csys_m_err,
            geN_ij.weight, geN_ij.npairs+geS_ij.npairs,
              ee_p_all, ee_m_all, np.sqrt(eeN.varxip+eeS.varxip), ge_p_all_ij, ge_m_all_ij, ge_p_all_ji, ge_m_all_ji, ede_p_all, dede_p_all, gde_p_all, sigma_gde_p_all,npairs_weighted ])
