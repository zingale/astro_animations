# compute the strength of Balmer lines in a star taking into account
# both excitation and ionization.  This follows the description from
# Carroll and Ostlie.

import math
import numpy
import pylab

def line_strength():

    ax = pylab.subplot(111)    
    pylab.subplots_adjust(left=0.175,right=0.95,bottom=0.15,top=0.95)

    # fundamental constants
    k_B    = 1.38e-16  # Boltzmann's constant [erg / K]
    m_e    = 9.11e-28  # electron mass [g]
    h      = 6.63e-27  # Planck's constant [erg s]
    eV2erg = 1.602e-12 # number of erg in 1 eV 
    xi_H   = 13.6      # ionization potential of H [eV]

    # electron pressure in the star (we assume that it is constant here)
    Pe = 200           # [dyn / cm^2]


    # define the temperature
    npts = 2000
    Tmax = 25000
    T = numpy.arange(npts, dtype=numpy.float64)*Tmax/(npts-1)


    # we are interested only in the ground state (n=1) and the first
    # excited state (n=2)
    
    # statistical weights in the H atom (g = 2n**2)
    g1 = 2
    g2 = 8
    
    # energy for n = 1 and n = 2 in the H atom
    E1 = -xi_H * eV2erg
    E2 = -xi_H/2**2 * eV2erg

    # partition functions for neutral H (ZI) and ionized H (ZII)
    # ZI = g1 + sum_{j=2}^{inf} g_j exp{-(E_j - E_1)/kT} 
    ZI = g1 + g2*numpy.exp(-(E2-E1)/(k_B*T))  # + ... (higher terms are negligible)
    ZII = 1                          # ionized H is just a proton


    #-------------------------------------------------------------------------

    # plot the fraction of electrons in the n = 2 state compared to all 
    # electrons (n = 1 and n = 2).  This is n_2 / (n_1 + n_2).  Here,
    # the ratio, n_2/n_1 is given by the Boltzmann equation.
    n2_over_n1 = (g2/g1)*numpy.exp(-(E2-E1)/(k_B*T))

    # f_ex = n_2 / (n_1 + n_2)
    f_ex = n2_over_n1 / (1.0 + n2_over_n1)

    pylab.plot(T, f_ex, 'r')

    pylab.xlabel(r"$T$ (K)")
    pylab.ylabel(r"$n_2/(n_1 + n_2)$")
    
    ax = pylab.gca()
    ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
    ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))

    ax.xaxis.major.formatter.set_powerlimits((-4, 5))
    ax.yaxis.major.formatter.set_powerlimits((-3, 4))

    f = pylab.gcf()
    f.set_size_inches(6.0,6.0)

    pylab.savefig("line_strength_excitation.png")


    #-------------------------------------------------------------------------
    
    # plot the fraction of H that is ionized compared to all the H (neutral
    # + ionized).  This is nII / (nI + nII).  Here, nII/nI is given by
    # the Saha equation.

    # following Carroll and Ostlie, we make the approximation that the 
    # electron pressure is constant.
    ne = Pe/(k_B*T)

    nII_over_nI = (2.0*ZII/ZI)*(1.0/ne) * \
        (2.0*math.pi*m_e*k_B*T/h**2)**1.5 * \
        numpy.exp(-xi_H*eV2erg/(k_B*T))

    # g_ion = nII / (nI + nII)
    g_ion = nII_over_nI / (1.0 + nII_over_nI)


    pylab.clf()

    pylab.plot(T, g_ion, 'r')

    pylab.xlabel(r"$T$ (K)")
    pylab.ylabel(r"$n_{II}/(n_{I} + n_{II})$")
    
    ax = pylab.gca()
    ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
    ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))

    ax.xaxis.major.formatter.set_powerlimits((-4, 5))
    ax.yaxis.major.formatter.set_powerlimits((-3, 4))

    f = pylab.gcf()
    f.set_size_inches(6.0,6.0)

    pylab.savefig("line_strength_ionization.png")


    #-------------------------------------------------------------------------
    
    # plot the fraction of H that is in level 2 compared to all the H (neutral
    # + ionized).  This is [n2 / (n1 + n2)] [nI / (nI + nII)].  

    h_tot = f_ex * (1.0 / (1.0 + nII_over_nI))

    pylab.clf()

    pylab.plot(T, h_tot, 'r')

    pylab.xlabel(r"$T$ (K)")
    pylab.ylabel(r"$n_{2}/n_{tot}$")
    
    ax = pylab.gca()
    ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
    ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))

    ax.xaxis.major.formatter.set_powerlimits((-4, 5))
    ax.yaxis.major.formatter.set_powerlimits((-3, 4))

    f = pylab.gcf()
    f.set_size_inches(6.0,6.0)

    pylab.savefig("line_strength.png")


if __name__== "__main__":
    line_strength()

