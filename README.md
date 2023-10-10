# CSK1000LF321

Cosmic shear analysis with updated shear measurement and calibration achieved using an updated version of the *lens*fit code and SKiLLS multi-band image simulations.

The repository documents the analysis code used in the publication: 

**KiDS-1000: Cosmology with improved cosmic shear measurements** ([Li et al. 2023](https://ui.adsabs.harvard.edu/abs/2023arXiv230611124L/abstract)).

It contains six parts:

- [A_Prepare_Shear_Cata](https://github.com/lshuns/CSK1000LF321/A_Prepare_Shear_Cata) for preparing shear and redshift catalogues from both KiDS data and SKiLLS simulations.

- [B_Shear_Bias](https://github.com/lshuns/CSK1000LF321/B_Shear_Bias) for calculating the mean shear biases for each tomographic bin.

- [C_Calc_2pt_Stats_Leiden](https://github.com/lshuns/CSK1000LF321/C_Calc_2pt_Stats_Leiden) for estimating the two-point shear correlation functions.

- [D_Covariance](https://github.com/lshuns/CSK1000LF321/D_Covariance) for computing the covariance matrix.

- [E_Cosmology](https://github.com/lshuns/CSK1000LF321/E_Cosmology) for performing cosmological inference.

- [paper_plot](https://github.com/lshuns/CSK1000LF321/paper_plot) for making the paper plots.