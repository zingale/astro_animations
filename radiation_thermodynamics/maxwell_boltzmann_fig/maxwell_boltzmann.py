# plot the Maxwell-Boltzmann distribution for H atoms at T = 10000 K

import math
import numpy
import pylab

def maxwell_boltzmann():

    ax = pylab.subplot(111)    

    k_B = 1.38e-16 # Boltzmann's constant [erg / K]
    m   = 1.67e-24 # mass of proton / H [g]
    T   = 10000    # temperature [K]

    npts = 2000
    vmax = 4.e6    # maximum velocity to plot [cm/s]

    v_mp = math.sqrt(2.0*k_B*T/m)
    v_rms = math.sqrt(3.0*k_B*T/m)

    v = numpy.arange(npts, dtype=numpy.float64)*vmax/(npts-1)

    MB = 4.0*math.pi*v**2 * \
        (m/(2.0*math.pi*k_B*T))**1.5 * \
        numpy.exp(-m*v**2/(2.0*k_B*T))

    pylab.plot(v, MB, 'r')

    pylab.plot([v_mp,v_mp], [0,7.e-7], 'k--')
    pylab.plot([v_rms,v_rms], [0,7.e-7 ], 'k--')

    pylab.xlabel("v (cm/s)")
    pylab.ylabel("P")
    
    ax = pylab.gca()
    ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
    ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))

    majorLocator = pylab.MultipleLocator(1.e6)
    ax.xaxis.set_major_locator(majorLocator)



    f = pylab.gcf()
    f.set_size_inches(6.0,6.0)

    pylab.savefig("maxwell_boltzmann.png")



if __name__== "__main__":
    maxwell_boltzmann()

