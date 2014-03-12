import math
import numpy
import pylab

def HR_mass():

    # solar temperature
    T_sun = 5777

    # temperature plotting range
    Tmin = 2000
    Tmax = 60000

    # luminosity plotting range
    Lmin = 1.e-4
    Lmax = 1.e6

    #-------------------------------------------------------------------
    # main sequence data taken from Carroll & Ostlie, Appendix G
    nstars = 6
    MS_M = numpy.zeros(nstars, numpy.float64)
    MS_T = numpy.zeros(nstars, numpy.float64)
    MS_R = numpy.zeros(nstars, numpy.float64)
    MS_L = numpy.zeros(nstars, numpy.float64)
    
    # spectral type
    MS_spectralTypes = ['E', 'F', 'H', 'I', 'J', 'M']
    MS_letter        = MS_spectralTypes
    # mass (solar masses)
    MS_M[:] = [17.5, 5.9, 2.9, 1,  0.51, 0.21]

    # temperature (K)
    MS_T[:] = [32500, 15200, 9800, T_sun, 3840, 3170]

    # radius (solar radii)
    MS_R[:] = [6.7, 3.2, 2.2, 1.0, 0.63, 0.29]

    # luminosity (solar luminosities)
    MS_L[:] = [32500, 480, 39.4, 1.0, 0.077, 0.0076]


    #-------------------------------------------------------------------
    # supergiants -- from C&O
    nstars = 4
    S_M = numpy.zeros(nstars, numpy.float64)
    S_T = numpy.zeros(nstars, numpy.float64)
    S_R = numpy.zeros(nstars, numpy.float64)
    S_L = numpy.zeros(nstars, numpy.float64)
    
    # spectral type
    S_spectralTypes = ['A', 'B', 'C', 'D']
    S_letter        = S_spectralTypes
    # mass (solar masses)
    S_M[:] = [25, 16, 10, 13]

    # temperature (K) -- T_sun sub to make things line up
    S_T[:] = [26000, 9730, T_sun, 3850]

    # radius (solar radii)
    S_R[:] = [25, 66, 190, 440]

    # luminosity (solar luminosities)
    S_L[:] = [260000, 35000, 30000, 38000]


    #-------------------------------------------------------------------
    # WDs -- made up
    nstars = 2
    W_T = numpy.zeros(nstars, numpy.float64)
    W_R = numpy.zeros(nstars, numpy.float64)
    W_L = numpy.zeros(nstars, numpy.float64)
    
    # spectral type
    W_letter        = ['K', 'L']

    # temperature (K)
    W_T[:] = [26000, 15000]

    # radius (solar radii)
    W_R[:] = [0.01, 0.01]

    # luminosity (solar luminosities)
    W_L[:] = [W_R[0]**2 * (W_T[0]/T_sun)**4, W_R[1]**2 * (W_T[1]/T_sun)**4]


    #-------------------------------------------------------------------
    # Sun, tip of giant branch -- made up
    nstars = 1
    G_T = numpy.zeros(nstars, numpy.float64)
    G_L = numpy.zeros(nstars, numpy.float64)
    
    # spectral type
    G_letter        = ['G']

    # temperature (K)
    G_T[:] = [4000]

    # luminosity (solar luminosities)
    G_L[:] = [300]


    ax = pylab.subplot(111)
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    # draw and label the main sequence
    pylab.scatter(MS_T,MS_L, s=100, marker="+", color="r", lw=2)
    #pylab.plot(T,L,'b-')
    n = 0
    while (n < len(MS_letter)):
        pylab.text(1.1*MS_T[n], MS_L[n], MS_letter[n], horizontalalignment="right")

        n += 1

    # draw and label the supergiants
    pylab.scatter(S_T,S_L, s=100, marker="+", color="r", lw=2)
    #pylab.plot(T,L,'b-')
    n = 0
    while (n < len(S_letter)):
        pylab.text(0.9*S_T[n], S_L[n], S_letter[n])

        n += 1


    # draw and label the WDs
    pylab.scatter(W_T,W_L, s=100, marker="+", color="r", lw=2)
    #pylab.plot(T,L,'b-')
    n = 0
    while (n < len(W_letter)):
        pylab.text(0.9*W_T[n], W_L[n], W_letter[n])

        n += 1


    # draw and label the giants
    pylab.scatter(G_T,G_L, s=100, marker="+", color="r", lw=2)
    n = 0
    while (n < len(G_letter)):
        pylab.text(0.9*G_T[n], G_L[n], G_letter[n])

        n += 1
 

    # draw x-axis (temperature) labels in a few spots
    Tplot = [40000, 20000, 10000, 5000, 2500]
    Tlabels = []
    n = 0
    while (n < len(Tplot)):
        Tlabels.append("%5.0f" % (Tplot[n]))
        n += 1
    
    locs, labels = pylab.xticks(Tplot, Tlabels)

    # reverse the x-axis
    pylab.axis([Tmin, Tmax, Lmin, Lmax])
    ax = pylab.gca()
    ax.set_xlim(ax.get_xlim()[::-1])

    # turn off minor ticks
    minorLocator = pylab.NullLocator()
    ax.xaxis.set_minor_locator(minorLocator)

    pylab.xlabel("$T$ (K)")
    pylab.ylabel(r"$L/L_\odot$")

    f = pylab.gcf()
    f.set_size_inches(6.0,6.0)

    pylab.savefig("HR_exam.eps")



if __name__== "__main__":
    HR_mass()

