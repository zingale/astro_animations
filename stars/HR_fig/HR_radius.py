# plot a H-R diagram.
# stellar properties taken from Carroll and Ostlie, Appendix G.
#
# optionally show lines of constant radius and the white dwarf
# population.
#
# M. Zingale (2010-12-13)


import numpy as np
import matplotlib.pyplot as plt

def HR_radius():

    # solar temperature
    T_sun = 5777.

    # temperature plotting range
    Tmin = 2000
    Tmax = 60000

    # luminosity plotting range
    Lmin = 1.e-4
    Lmax = 1.e6

    # draw radius lines? (1 = yes, 0 = no)
    draw_radii = 1

    # draw the WD region?
    draw_WDs = 1

    # label the stars by mass?
    label_masses = 1

    #-------------------------------------------------------------------------
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

    #-------------------------------------------------------------------------


    ax = plt.subplot(111)
    ax.set_xscale('log')
    ax.set_yscale('log')

    # draw and label the main sequence
    plt.scatter(T, L, s=100, marker="+", color="r", lw=2)
    plt.plot(T, L, 'b-')
    for mstar, tstar, lstar, spec_star in zip(M, T, L, spectral_types):

        if label_masses:
            if mstar >= 1.0:
                plt.text(1.05*tstar, lstar, fr"{mstar:4.1f} $M_\odot$", fontsize=10,
                         horizontalalignment="right", verticalalignment="center", color="r")
            else:
                plt.text(1.05*tstar, lstar, fr"{mstar:3.2f} $M_\odot$", fontsize=10,
                         horizontalalignment="right", verticalalignment="center", color="r")

        plt.text(0.9*tstar, lstar, spec_star)


    if draw_radii:

        # draw in lines of constant radius
        Rvals = [0.0001, 0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000, 10000]

        nTpts = 25
        dlogT = (np.log10(Tmax) - np.log10(Tmin))/(nTpts-1)

        Tplot = 10.0**(np.arange(nTpts)*dlogT + np.log10(Tmin))

        for n in range(len(Rvals)):

            # L/L_sun = (R/R_sun)**2 (T/T_sun)**4
            # compute L in L_sun units
            L = Rvals[n]**2 * (Tplot/T_sun)**4

            plt.plot(Tplot, L, linestyle='--', color="0.5")

            # draw in labels -- here we space the T position of the labels equally
            # in log, by empericallly picking T = 40000 for L = 1.e-4 and
            # T = 3000 for L = 1.e4, and doing a line in log-space
            slope = (np.log10(40000) - np.log10(3000))/(np.log10(1.e-4) - np.log10(1.e4))
            xpos = np.log10(40000) + slope*(np.log10(Rvals[n]) - np.log10(1.e-4))
            xpos = 10.0**xpos

            ypos = Rvals[n]**2 * (xpos/T_sun)**4


            if (xpos > Tmin and xpos < Tmax and ypos > Lmin and ypos < Lmax):
                plt.text(xpos, ypos, fr"{Rvals[n]:g} $R_\odot$", color="0.5")

    if draw_WDs:

        # draw a box along the R = 0.01 R_sun line to indicate the approximate
        # position of the white dwarfs
        L3 = 0.01**2 * (30000./T_sun)**4
        L75 = 0.01**2 * (7500./T_sun)**4
        plt.fill([30000, 30000, 7500, 7500], [1.5*L3, 0.66*L3, 0.66*L75, 1.5*L75],
                 alpha=0.20, facecolor="b")

        L15 = 0.01**2 * (15000./T_sun)**4
        plt.text(15000, L15, "white dwarfs", color="b", horizontalalignment="center")

    # draw x-axis (temperature) labels in a few spots
    Tplot = [40000, 20000, 10000, 5000, 2500]
    Tlabels = []

    for n in range(len(Tplot)):
        Tlabels.append("%5.0f" % (Tplot[n]))

    locs, labels = plt.xticks(Tplot, Tlabels)

    # reverse the x-axis
    plt.axis([Tmin, Tmax, Lmin, Lmax])
    ax = plt.gca()
    ax.set_xlim(ax.get_xlim()[::-1])

    # turn off minor ticks
    minorLocator = plt.NullLocator()
    ax.xaxis.set_minor_locator(minorLocator)

    plt.xlabel("$T$ (K)")
    plt.ylabel(r"$L/L_\odot$")

    plt.tight_layout()

    f = plt.gcf()
    f.set_size_inches(6.0, 7.0)

    plt.savefig("HR_radius.png")


if __name__ == "__main__":
    HR_radius()
