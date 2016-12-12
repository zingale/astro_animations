#!/bin/env python

# plot the periods vs. semi-major axes for the planets in our solar
# system on a log-log plot, and draw a line through it showing that
# that obey Kepler's third law.
#
# optionally plot Pluto and the Galilean moons (the latter shows that
# the constant in Kepler's third law is different for bodies orbiting
# Jupiter.
#
# M. Zingale (2010-12-13)

import numpy as np
import matplotlib.pyplot as plt

def harmonic_law():

    plot_moons = False
    plot_pluto = False
    
    plot_extend_axes = True

    names = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", 
             "Uranus", "Neptune", "Pluto"]

    names_J = ["Io", "Europa", "Ganymede", "Callisto"]

    a = np.zeros(9, np.float64)
    P = np.zeros(9, np.float64)

    a_J = np.zeros(4, np.float64)
    P_J = np.zeros(4, np.float64)

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
    if plot_moons or plot_extend_axes:
        min_a = 0.001
        min_P = 0.001
    else:
        min_a = 0.1
        min_P = 0.1

    

    ax = plt.subplot(111)
    ax.set_xscale('log')
    ax.set_yscale('log')


    # solar system
    if plot_pluto:
        plt.scatter(a, P, s=100, marker="+", color="r", lw=2)
    else:
        plt.scatter(a[0:8], P[0:8], s=100, marker="+", color="r", lw=2)

    atemp = np.arange(2)*(2*np.max(a) - 0.5*np.min(a)) + 0.5*np.min(a)
    plt.plot(atemp, atemp**1.5, 'k--')

    if plot_pluto:
        num = 9
    else:
        num = 8

    for i in range(num):
        plt.text(a[i]*0.8,P[i], names[i],horizontalalignment="right", 
                   color="r", fontsize=10)


    if plot_moons:

        # Galilean moons
        plt.scatter(a_J,P_J, s=50, marker="o", color="b", lw=2)

        atemp = np.arange(2)*(2*np.max(a_J) - 0.5*np.min(a_J)) + 0.5*np.min(a_J)
        plt.plot(atemp, P_J[0]*(atemp/a_J[0])**1.5, 'k--')

        for i in range(4):
            plt.text(a_J[i]*1.5,P_J[i], names_J[i],horizontalalignment="left",
                       verticalalignment="center",
                       color="b",fontsize=10)


    plt.xlabel("semi-major axis (AU)")
    plt.ylabel("period (yr)")
        
    plt.axis([min_a,100,min_P,1000])

    plt.title(r"Kepler's Third Law",fontsize=11)

    f = plt.gcf()
    f.set_size_inches(7.0,6.0)

    plt.tight_layout()

    outbase = "harmonic_law"
    if plot_pluto: outbase += "_pluto"
    if plot_moons: outbase += "_moons"
        
    plt.savefig(outbase + ".png")
    plt.savefig(outbase + ".eps")



if __name__== "__main__":
    harmonic_law()

