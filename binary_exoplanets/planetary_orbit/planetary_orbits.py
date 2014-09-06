import math
import numpy
import matplotlib.pyplot as plt

import anim_solvers.binary_integrator as bi

# compute the orbits of planet around a star -- this is based off of
# binary_stars.py
#
# Here we put the center of mass at the origin.  
#
# This version allows for elliptical orbits with some arbitrary
# orientation wrt to the observer (although, still face-on)
#
# The plotting assumes that M_2 < M_1
#  
# M. Zingale (2011-09-12)

# we work in CGS units
G = bi.G
M_sun = bi.M_sun
AU = bi.AU
year = bi.year


def doit():

    # set the masses
    M_star1 = 50*M_sun           # star 1's mass
    M_star2 = M_sun      # star 2's mass

    # set the semi-major axis of the star 2 (and derive that of star 1)
    # M_star2 a_star2 = -M_star1 a_star1 (center of mass)
    a_star2 = 1.0*AU
    a_star1 = (M_star2/M_star1)*a_star2  

    # set the eccentricity
    ecc = 0.0

    # set the angle to rotate the semi-major axis wrt the observer
    theta = math.pi/6.0

    # create the solar system container
    b = bi.Binary(M_star1, M_star2, a_star1+a_star2, ecc, theta)

    
    # set the timestep in terms of the orbital period
    dt = b.P/360.0        
    tmax = 2.0*b.P  # maximum integration time

    s1, s2 = b.integrate(dt, tmax)


    # ================================================================
    # plotting
    # ================================================================

    iframe = 0

    for n in range(len(s1.t)):

        plt.clf()

        plt.subplots_adjust(left=0.1,right=0.9,bottom=0.1,top=0.9)

        a = plt.gca()
        a.set_aspect("equal", "datalim")
        plt.axis("off")

        # if e = 0 and M_star1 = M_star2, then the orbits lie on top of one
        # another, so plot only a single orbital line.

        # plot star 1's orbit and position
        symsize = 500
        if not (b.M1 == b.M2 and b.e == 0.0):
            plt.plot(s1.x, s1.y, color="r")
        else:
            plt.plot(s1.x, s1.y, color="k")

        plt.scatter([s1.x[n]], [s1.y[n]], s=symsize, color="r")

        # plot star 2's orbit and position
        symsize = 1000*(b.M2/b.M1)
        if not (b.M1 == b.M2 and b.e == 0.0):
            plt.plot(s2.x, s2.y, color="g")

        plt.scatter([s2.x[n]], [s2.y[n]], s=symsize, color="g")

        plt.scatter([0], [0], s=150, marker="x",color="k", zorder=1000)

        plt.axis([-1.4*a_star2, 1.4*a_star2, -1.4*a_star2, 1.4*a_star2])

        f = plt.gcf()
        f.set_size_inches(7.2,7.2)
  
        plt.savefig("planetary_orbit_%04d.png" % iframe)

        iframe += 1

    
if __name__== "__main__":
    doit()


    
        
