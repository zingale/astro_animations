#!/bin/env python

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, NullLocator

# draw blackbody curves of varying temperature
# here we plot I_lambda

# M. Zingale 

# physical constants (cgs)
h = 6.626068e-27
c = 2.99792458e10
k = 1.3806503e-16


def I_lambda(l, T):
    I = (2.0*h*c*c/(l**5))/(np.exp(h*c/(l*k*T)) - 1.0)
    return I


def blackbody():

    # define wavelength
    npts = 250
    l = np.logspace(-14.0, 5.0, npts, endpoint=True)

    npts_T = 500
    T = np.logspace(0.0, 10.0, npts_T, endpoint=True)


    # some temperatures to highlight
    T_ref = np.array([1.e2, 1.e4, 1.e6, 1.e8])
    I_ref = []

    for T0 in T_ref:
        I_ref.append(I_lambda(l, T0))

    I_ref = np.array(I_ref)


    for n in range(len(T)):
        print n, T[n]

        I = I_lambda(l, T[n])


        # plotting
        plt.clf()

        plt.subplots_adjust(left=0.125,right=0.75,
                              bottom=0.1,top=0.9,wspace=0.1)       
        
        plt.loglog(l,I,color="g",linewidth=2.0)

        for T0, I0 in zip(T_ref, I_ref):

            # reference lines
            if T[n] >= T0:

                plt.loglog(l,I0, color="#888888", linewidth=2.0, alpha=0.5)

                # use Wien's law to find the max
                lmax = 0.29/T0
                Imax = I_lambda(lmax, T0)

                exp = int(math.log10(T0))
                pre = T0/10.0**exp
                if pre == 1.0:
                    plt.text(1.1*lmax,1.1*Imax, 
                             r"$T = 10^{%d} \, \mathrm{K} $" % (exp), 
                             color="#666666")
                else:
                    plt.text(1.1*lmax,1.1*Imax, 
                             r"$T = %3.2f \times 10^{%d} \, \mathrm{K} $" % (pre, exp), 
                             color="#666666")



        # visible wavelengths: 400 nm to 700 nm  (4e-5 cm to 7e-5 cm)
        Imin = 1.e-6
        Imax = 1.e48

        red = 7.e-5
        blue = 4.e-5

        plt.fill([blue, blue, red,  red,  blue],
                 [Imin, Imax, Imax, Imin, Imin], alpha=0.20, facecolor="b")

        plt.axis([np.min(l), np.max(l), Imin, Imax])

        plt.xlabel("wavelength [cm]")
        plt.ylabel("Intensity")
        plt.title("Blackbody Radiation")



        # draw a thermometer 
        plt.axes([0.885,0.1,0.015,0.8])

        # we do a lot of thin lines here, instead of using a thick
        # line width, because a thick line extends vertically as
        # well as horizontally, and does not line up to the right place
        # on the axis
        therm_xmax = 0.01
        iter = 0
        while (iter < 10):
            plt.semilogy([iter*therm_xmax/10.0,iter*therm_xmax/10.0],
                           [T[0],T[n]], color='r',linewidth=2)
            iter += 1


        plt.axis([0,therm_xmax,T[0],T[npts_T-1]])

        axis = plt.gca()
        #majorLocator = MultipleLocator(20)
        majorLocator = NullLocator()
        axis.xaxis.set_major_locator(majorLocator)

        plt.ylabel("temperature [K]")

        fig = plt.gcf()
        fig.set_size_inches(9.6,7.2)

        plt.savefig("blackbody_%03d.png" % n)


if __name__== "__main__":
    blackbody()


    
        
