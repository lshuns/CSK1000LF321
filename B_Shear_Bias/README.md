# Calculate Shear measurement biases

It heavily depends on the [biasEstimation](https://github.com/KiDS-WL/MultiBand_ImSim/tree/main/biasEstimation) and [correction_varShear_PSFmodelling](https://github.com/KiDS-WL/MultiBand_ImSim/tree/main/correction_varShear_PSFmodelling) repositories.

[A_run_raw_m_ZB_noRewei.sh](https://github.com/lshuns/CSK1000LF321/tree/main/B_Shear_Bias/A_run_raw_m_ZB_noRewei.sh)

- Calculate the raw shear biases from the constant shear image simulations without data re-weighting.

[A_run_raw_m_ZB.sh](https://github.com/lshuns/CSK1000LF321/tree/main/B_Shear_Bias/A_run_raw_m_ZB.sh)

- Calculate the raw shear biases from the constant shear image simulations with data re-weighting.

[B_blending_fraction_2D.py](https://github.com/lshuns/CSK1000LF321/tree/main/B_Shear_Bias/B_blending_fraction_2D.py)

- Estimate the blending fraction. This is derived from [A_blending_fraction_2D.py](https://github.com/KiDS-WL/MultiBand_ImSim/blob/main/correction_varShear_PSFmodelling/A_blending_fraction_2D.py).

[C_dm_inBlendingSample_2D.py](https://github.com/lshuns/CSK1000LF321/tree/main/B_Shear_Bias/C_dm_inBlendingSample_2D.py)

- Calculate the shear bias difference between constant shear simulations and variable shear simulations using the blending-only samples. This is derived from [B_dm_inBlendingSample_2D.py](https://github.com/KiDS-WL/MultiBand_ImSim/blob/main/correction_varShear_PSFmodelling/B_dm_inBlendingSample_2D.py).

[D_apply_dm_varShear_2D.py](https://github.com/lshuns/CSK1000LF321/tree/main/B_Shear_Bias/D_apply_dm_varShear_2D.py)

- Adjust the raw shear biases to account for the dm estimated from variable shear simulations. This is derived from [C_apply_dm_varShear_2D.py](https://github.com/KiDS-WL/MultiBand_ImSim/blob/main/correction_varShear_PSFmodelling/C_apply_dm_varShear_2D.py).

[E_dm_PSFmodelling.py](https://github.com/lshuns/CSK1000LF321/tree/main/B_Shear_Bias/E_dm_PSFmodelling.py)

- Calculate the shear bias difference between idealised simulations and simulations with PSF modeling procedures. This is derived from [D_dm_PSFmodelling.py](https://github.com/KiDS-WL/MultiBand_ImSim/blob/main/correction_varShear_PSFmodelling/D_dm_PSFmodelling.py).

[F_apply_dm_PSFmodelling.py](https://github.com/lshuns/CSK1000LF321/tree/main/B_Shear_Bias/F_apply_dm_PSFmodelling.py)

- Correct the shear biases to account for the dm estimated from the PSF modeling test. This is derived from [E_apply_dm_PSFmodelling.py](https://github.com/KiDS-WL/MultiBand_ImSim/blob/main/correction_varShear_PSFmodelling/E_apply_dm_PSFmodelling.py).

[Z_alpha_c.py](https://github.com/lshuns/CSK1000LF321/tree/main/B_Shear_Bias/Z_alpha_c.py)

- Calculate the PSF leakage.