# illustrate the Gamov peak as well as an approximate Gaussian fit to it.
# compare to Clayton Fig. 4-7

import numpy as np
import matplotlib.pyplot as plt

# this is all based on Clayton, Ch. 4.

# nuclei
Z_1 = 1.0
A_1 = 1.0
Z_2 = 6.0
A_2 = 12.0

A = A_1*A_2/(A_1 + A_2)

T = 3.e7   # K

# constants
k = 1.38e-16   # erg/K
keV_erg = 1.602e-9  # erg/keV
k_keV = k/keV_erg   # keV/K

b = 31.2*np.sqrt(A)*Z_1*Z_2   # Clayton Eq. 4-41 units: keV**1/2
E_0 = (b*k_keV*T/2.0)**(2./3.)  # Clayton Eq. 4-46


def gamow(E, T):
    # E should be in keV
    return np.exp(-E/(k_keV*T) - b/np.sqrt(E))


def max_boltz(E, T):
    return np.exp(-E/(k_keV*T))


def tunnel(E, T):
    return np.exp(- b/np.sqrt(E))


def gaussian_fit(E, T):
    # this is Clayton Eqs. 4-48 to 4-51
    tau = 3.0*E_0/(k_keV*T)
    C = np.exp(-tau)
    Delta = 4.0/np.sqrt(3) * np.sqrt(E_0 * k_keV*T)

    return C*np.exp(-((E-E_0)*2/Delta)**2)


fig, ax = plt.subplots()

E_0 = (b*k_keV*T/2.0)**(2./3.)

print("E_0 = ", E_0)

E = np.linspace(0, 100, 200)

ax.plot(E, 1.e15*gamow(E, T),
        label=r"$e^{-E/kT - bE^{-1/2}} (\times 10^{15})$")
ax.plot(E, max_boltz(E, T),
        label=r"Maxwell-Boltzmann tail")
ax.plot(E, 1.e7*tunnel(E, T),
        label=r"tunneling factor $(\times 10^7)$")

ax.plot([E_0, E_0], [0, 1], ls=":", color="0.5")

ax.legend(frameon=False, fontsize=12)

ax.set_xlim(20, 60)
ax.set_ylim(0, 1.e-4)
ax.yaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))


ax.set_xlabel("E (keV)")

ax.set_title(r"Gamow peak: $^{12}\mathrm{C}(p,\gamma){^{13}\mathrm{N}},\ T = 3\times 10^7~\mathrm{K}$", fontsize=12)

fig.tight_layout()
fig.set_size_inches(7.2, 7.2)

plt.savefig("gamow.png")
