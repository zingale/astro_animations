import math
import numpy as np
import matplotlib.pyplot as plt

import anim_solvers.binary_integrator as bi

# compute the orbits of stars in a binary system.
#
# Here we put the center of mass at the origin.  
#
# This version allows for elliptical orbits with some arbitrary
# orientation wrt to the observer (although, still face-on)
#
# The plotting assumes that M_2 < M_1
#  
# M. Zingale (2009-02-12)

# we work in CGS units
G = bi.G
M_sun = bi.M_sun
AU = bi.AU
year = bi.year


def doit():

    # set the masses
    M_star1 = M_sun           # star 1's mass
    M_star2 = 0.6*M_sun      # star 2's mass

    # set the semi-major axis of the star 2 (and derive that of star 1)
    # M_star2 a_star2 = -M_star1 a_star1 (center of mass)
    a_star2 = 1.0*AU
    a_star1 = (M_star2/M_star1)*a_star2  

    # set the eccentricity
    ecc = 0.4

    # set the angle to rotate the semi-major axis wrt the observer
    theta = math.pi/6.0

    # create the solar system container
    b = bi.Binary(M_star1, M_star2, a_star1 + a_star2, ecc, theta)

    
    # set the timestep in terms of the orbital period
    dt = b.P/360.0        
    tmax = 2.0*b.P  # maximum integration time

    s1, s2 = b.integrate(dt, tmax)


    # ================================================================
    # plotting
    # ================================================================

    plt.clf()

    plt.subplots_adjust(left=0.1,right=0.9,bottom=0.1,top=0.9)

    a = plt.gca()
    a.set_aspect("equal", "datalim")
    plt.axis("off")

    n = 0
    
    plt.scatter([0], [0], s=150, marker="x", color="k")

    # plot star 1's orbit and position
    symsize = 200
    plt.plot(s1.x, s1.y, color="r")
    

    # plot star 2's orbit and position
    #symsize = 200*(M_star2/M_star1)
    plt.plot(s2.x, s2.y, color="g", linestyle="--")

    plt.scatter([s2.x[n]], [s2.y[n]], s=symsize, color="g", marker='h')

    # reference points -- in terms of total number of steps integrated
    f1 = 0.0943
    f2 = 0.25
    f3 = 0.5 - f1

    npts = len(s2.t)

    xc = 0.5*(s2.x[0] + s2.x[0.25*npts])
    yc = 0.5*(s2.y[0] + s2.y[0.25*npts])
    
    # compute the angle of f1 wrt the initial position (just for debugging)
    a1 = math.atan2(s2.y[f1*npts] - yc, s2.x[f1*npts] - xc)
    a0 = math.atan2(s2.y[0] - yc, s2.x[0] - xc)

    print (a1-a0)*180./math.pi

    plt.scatter([s2.x[f1*npts]], [s2.y[f1*npts]], 
                s=symsize, color="g", marker='h')
    plt.scatter([s2.x[f2*npts]], [s2.y[f2*npts]], 
                s=symsize, color="g", marker='h')
    plt.scatter([s2.x[f3*npts]], [s2.y[f3*npts]], 
                s=symsize, color="g", marker='h')

    # label the points
    plt.text(s2.x[n]*1.15, s2.y[n], "A4")
    plt.text(s2.x[f1*npts]*1.15, s2.y[f1*npts]*1.15, "A1")
    plt.text(s2.x[f2*npts]*1.15, s2.y[f2*npts], "A2")
    plt.text(s2.x[f3*npts]*0.85, s2.y[f3*npts]*1.15, "A3")


    plt.axis([-1.4*b.a2,0.8*b.a2,-1.4*b.a2,0.8*b.a2])

    f = plt.gcf()
    f.set_size_inches(7.2,7.2)
  
    plt.savefig("binary_fig.eps", bbox_inches="tight")

    
if __name__== "__main__":
    doit()


    
        
