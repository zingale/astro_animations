# plot the periods vs. semi-major axes for the planets in our solar
# system on a log-log plot, and draw a line through it showing that
# that obey Kepler's third law.
#
# optionally plot Pluto and the Galilean moons (the latter shows that
# the constant in Kepler's third law is different for bodies orbiting
# Jupiter.
#
# M. Zingale (2010-12-13)

import math
import numpy
import pylab

def harmonic_law():

    plot_moons = 1
    plot_pluto = 1

    names = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]

    names_J = ["Io", "Europa", "Ganymede", "Callisto"]

    a = numpy.zeros(9, numpy.float64)
    P = numpy.zeros(9, numpy.float64)

    a_J = numpy.zeros(4, numpy.float64)
    P_J = numpy.zeros(4, numpy.float64)

    # solar system
    # semi-major axes
    a[:] = [0.39, 0.72, 1.00, 1.52,  5.20,  9.54, 19.22,  30.06, 39.48]

    # period
    P[:] = [0.24, 0.62, 1.00, 1.88, 11.86, 29.46, 84.01, 164.8,  248.09]


    # Jovian moons
    # semi-major axes (km)
    a_J[:] = [4.218e5, 6.711e5, 1.0704e6, 1.8827e6]
    a_J[:] = a_J[:]/1.49598e8 # (km / AU conversion)

    # period (days)
    P_J[:] = [1.769, 3.551, 7.155, 16.69]
    P_J[:] = P_J[:]/365.24 # (days / yr)

    
    # minimum semi-major axis value on x-axis
    if (plot_moons == 1):
        min_a = 0.001
        min_P = 0.001
    else:
        min_a = 0.1
        min_P = 0.1

    

    ax = pylab.subplot(111)
    ax.set_xscale('log')
    ax.set_yscale('log')


    # solar system
    pylab.scatter(a,P, s=100, marker="+", color="r", lw=2)

    atemp = numpy.arange(2)*(2*numpy.max(a) - 0.5*numpy.min(a)) + 0.5*numpy.min(a)
    pylab.plot(atemp, atemp**1.5, 'k--')

    i = 0
    if (plot_pluto == 1):
        num = 9
    else:
        num = 8

    while (i < num):
        pylab.text(a[i]*0.8,P[i], names[i],horizontalalignment="right", 
                   color="r", fontsize=10)
        i += 1


    print a_J
    print P_J

    if (plot_moons == 1):

        # Galilean moons
        pylab.scatter(a_J,P_J, s=50, marker="o", color="b", lw=2)

        atemp = numpy.arange(2)*(2*numpy.max(a_J) - 0.5*numpy.min(a_J)) + 0.5*numpy.min(a_J)
        pylab.plot(atemp, P_J[0]*(atemp/a_J[0])**1.5, 'k--')

        i = 0
        while (i < 4):
            pylab.text(a_J[i]*1.5,P_J[i], names_J[i],horizontalalignment="left",
                       verticalalignment="center",
                       color="b",fontsize=10)
            i += 1


    pylab.xlabel("semi-major axis (AU)")
    pylab.ylabel("period (yr)")

        
    pylab.axis([min_a,100,min_P,1000])

    pylab.title(r"Kepler's Third Law",fontsize=11)

    f = pylab.gcf()
    f.set_size_inches(7.0,6.0)

    pylab.savefig("harmonic_law.png")
    pylab.savefig("harmonic_law.eps")



if __name__== "__main__":
    harmonic_law()

