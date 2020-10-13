import numpy as np
import matplotlib.pyplot as plt

import stellar_properties as sp


fig = plt.figure()
ax = fig.add_subplot(111)

plt.scatter(sp.M, sp.R, s=100, marker="+", color="C0", lw=2)
ax.plot(sp.M, sp.R, color="C0")

for mstar, rstar, spec_star in zip(sp.M, sp.R, sp.spectral_types):

    if mstar >= 1.0:
        plt.text(0.95*mstar, rstar, fr"{mstar:4.1f} $M_\odot$", fontsize=10,
                 horizontalalignment="right", verticalalignment="center", color="C0")
    else:
        plt.text(0.95*mstar, rstar, fr"{mstar:3.2f} $M_\odot$", fontsize=10,
                 horizontalalignment="right", verticalalignment="center", color="C0")

    plt.text(1.2*mstar, rstar, spec_star)

ax.set_ylabel(r"$R/R_\odot$")
ax.set_xlabel(r"$M/M_\odot$")
ax.set_xscale("log")
ax.set_yscale("log")

# theoretical scalings

# high mass, R ~ M**15/19
x = np.logspace(np.log10(1.0), np.log10(20.0), 10, endpoint=True)
ax.plot(x, x**(15./19), label=r"$R \sim M^{3/7}$", color="C1", linestyle=":")

ax.legend(frameon=False)

fig.savefig("ms_mass_vs_radius.png")
