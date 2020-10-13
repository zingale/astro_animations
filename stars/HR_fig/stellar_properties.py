import numpy as np

# solar temperature
T_sun = 5777.

# main sequence data taken from Carroll & Ostlie, Appendix G
nstars = 11
M = np.zeros(nstars, np.float64)
T = np.zeros(nstars, np.float64)
R = np.zeros(nstars, np.float64)
L = np.zeros(nstars, np.float64)

# spectral type
spectral_types = ['O8', 'B0', 'B3', 'B5', 'A0', 'A5', 'F5', 'Sun', 'K5', 'M0', 'M5']

# mass (solar masses)
M[:] = [23, 17.5, 7.6, 5.9, 2.9, 1.8, 1.2, 1, 0.67, 0.51, 0.21]

# temperature (K)
T[:] = [35800, 32500, 18800, 15200, 9800, 8190, 6650, T_sun, 4410, 3840, 3170]

# radius (solar radii)
R[:] = [10.0, 6.7, 3.8, 3.2, 2.2, 1.8, 1.2, 1.0, 0.80, 0.63, 0.29]

# luminosity (solar luminosities)
L[:] = [147000, 32500, 1580, 480, 39.4, 12.3, 2.56, 1.0, 0.216, 0.077, 0.0076]



