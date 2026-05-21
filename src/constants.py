import numpy as np

GAMMA_FIXED = 9.6e-8
K_PERM_M2 = 2.39e-13
MU_VISC_PAS = 0.7e-3
K_DARCY_FIXED = (K_PERM_M2 * 100.0) / (MU_VISC_PAS / 3600.0)
H_D_FIXED = 1400.0
INV_L_DURA = 0.01  # L=1mm char length
phi_CSF = 0.7
D_D = 3.5e-5

TIME_POINTS = np.array([0.0, 4.0, 24.0, 48.0, 72.0])
