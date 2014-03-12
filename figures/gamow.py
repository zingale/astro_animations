# illustrate the Gamov peak as well as an approximate Gaussian fit to it.
# compare to Clayton Fig. 4-7

import numpy
import math
import pylab

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

b = 31.2*math.sqrt(A)*Z_1*Z_2   # Clayton Eq. 4-41 units: keV**1/2
E_0 = (b*k_keV*T/2.0)**(2./3.)  # Clayton Eq. 4-46 

def gamow(E, T):
    # E should be in keV
    return numpy.exp(-E/(k_keV*T) - b/numpy.sqrt(E))


def gaussian_fit(E, T):
    # this is Clayton Eqs. 4-48 to 4-51
    tau = 3.0*E_0/(k_keV*T)
    C = math.exp(-tau)
    Delta = 4.0/math.sqrt(3) * math.sqrt(E_0 * k_keV*T)

    return C*numpy.exp(-((E-E_0)*2/Delta)**2)


E_0 = (b*k_keV*T/2.0)**(2./3.)

print "E_0 = ", E_0

E = numpy.linspace(0, 100, 200)


pylab.plot(E, gamow(E,T), label=r"$e^{-E/kT - bE^{-1/2}}$")
pylab.plot(E, gaussian_fit(E,T), color="r", ls=":", label=r"$e^{-\tau} \cdot e^{-[(E-E_0)/(\Delta/2)]^2}$")

pylab.legend(frameon=False, fontsize=12)

ax = pylab.gca()
ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))


pylab.xlabel("E (keV)")
pylab.title(r"Gamow peak + approximate Gaussian: $^{12}\mathrm{C}(p,\gamma){^{13}\mathrm{N}},\ T = 3\times 10^7~\mathrm{K}$", fontsize=12)

pylab.tight_layout()

pylab.savefig("gamow.png")
