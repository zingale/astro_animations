"""
a simple figure showing a planetary transit
"""

import numpy as np
import matplotlib.pyplot as plt

# we work in CGS units
G = 6.67428e-8        # cm^3 g^{-1} s^{-2}
M_sun = 1.98892e33    # g
R_sun = 7.e10         # cm
L_sun = 4.e33         # erg/s
AU = 1.49598e13       # cm
year = 3.1556926e7    # s

# planet properties
a = 0.5*AU
R_p = 0.0175*R_sun

# star properties (K0 star)
M_star = 0.78*M_sun
R_star = 0.85*R_sun
L_star = 0.4*L_sun

# find the period
P = np.sqrt(4.0*np.pi**2*a**3/(G*M_star))

print(f"period = {P/year} yr")

# find the velocity of the planet in the orbit
v = 2.0*np.pi*a/P

# dimming during transit
df = (R_p/R_star)**2

f_normal = 1.0
f_transit = f_normal - df

print("fluxes = ", f_normal, f_transit)

# length of transit
t_transit = R_star/v

# let's model 2 t_transit, with the transit centered in this interval
t_min = 0.0
t_max = 2.0*t_transit

# when does the transit begin
t_mid = 0.5*(t_max - t_min)
t_begin = t_mid - t_transit/2
t_end = t_mid + t_transit/2

# smooth by the radius crossing of the planet
dt_planet = R_p/v

t_begin_a = t_begin - 0.5*dt_planet
t_begin_b = t_begin + 0.5*dt_planet

t_end_a = t_end - 0.5*dt_planet
t_end_b = t_end + 0.5*dt_planet

# now draw it!
npts = 1000

t = np.linspace(t_min, t_max, npts)

f = f_normal*np.ones_like(t)

ib = np.logical_and(t > t_begin_a, t < t_begin_b)
f[ib] = f_normal + (f_transit - f_normal) / \
    (t_begin_b - t_begin_a)*(t[ib]-t_begin_a)

f[np.logical_and(t > t_begin_b, t < t_end_a)] = f_transit

ie = np.logical_and(t > t_end_a, t < t_end_b)
f[ie] = f_transit + (f_normal - f_transit)/(t_end_b - t_end_a)*(t[ie]-t_end_a)

fig, ax = plt.subplots()

ax.plot((t-t_mid)/3600.0, f)

ax.set_ylim(0.9995, 1.0001)

# turn off the offset text on the y-axis
ax.get_yaxis().get_major_formatter().set_useOffset(False)

ax.set_xlabel(r"$t$ (hours)", fontsize="large")
ax.set_ylabel(r"$f/f_\star$", fontsize="large")

fig.savefig("transit.png")
fig.savefig("transit.pdf")

