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


def astrometric_binary():

    # set the masses
    M_star1 = M_sun           # star 1's mass
    M_star2 = 0.5*M_sun      # star 2's mass

    # set the semi-major axis of the star 2 (and derive that of star 1)
    # M_star2 a_star2 = -M_star1 a_star1 (center of mass)
    a_star2 = 1.0*AU
    a_star1 = (M_star2/M_star1)*a_star2  

    # set the eccentricity
    ecc = 0.0

    # set the angle to rotate the semi-major axis wrt the observer
    theta = math.pi/2.0

    b = bi.Binary(M_star1, M_star2, a_star1+a_star2, ecc, theta)

    # velocity of the center of mass
    v = 10*a_star2/b.P

    
    # set the timestep in terms of the orbital period
    dt = b.P/360.0        
    tmax = 2.0*b.P  # maximum integration time

    s1, s2 = b.integrate(dt, tmax)


    # now move the center of mass according to the velocity defined
    # above, and update the positions
    x_cm = np.zeros(len(s1.t), np.float64)
    y_cm = np.zeros(len(s2.t), np.float64)

    for n in range(len(s1.t)):

        x_cm[n] = v*s1.t[n]
        y_cm[n] = 0.0

        s1.x[n] += v*s1.t[n]
        s2.x[n] += v*s2.t[n]


    # ================================================================
    # plotting
    # ================================================================

    plt.clf()

    plt.subplots_adjust(left=0.1,right=0.9,bottom=0.1,top=0.9)

    a = plt.gca()
    a.set_aspect("equal", "datalim")
    plt.axis("off")

    # plot the center of mass motion
    plt.scatter([0], [0], s=50, marker="x", color="k")
    plt.plot(x_cm, y_cm, color="0.5", linestyle="--")

    # plot star 1's orbit and position
    plt.plot(s1.x, s1.y, color="r")
    plt.scatter([s1.x[0]],[s1.y[0]], s=200, color="r")

    # plot star 2's orbit and position
    plt.plot(s2.x, s2.y, color="g")
    plt.scatter([s2.x[0]], [s2.y[0]], s=100, color="g")
        
    plt.text(10*a_star2, 3*a_star2, "mass ratio: %3.2f" % (b.M1/b.M2), 
             color="k", verticalalignment="center")
    plt.text(10*a_star2, 2.5*a_star2, "eccentricity: %3.2f" % (ecc), 
             color="k", verticalalignment="center")

    # labels
    plt.plot([0.1*a_star2,1.1*a_star2], [2.9*a_star2,2.9*a_star2], 
             color="0.5", linestyle="--")
    plt.text(1.3*a_star2, 2.9*a_star2, "path of center of mass",
             verticalalignment="center", color="0.5")

    plt.plot([0.1*a_star2,1.1*a_star2], [2.4*a_star2,2.4*a_star2], color="r")
    plt.text(1.3*a_star2, 2.4*a_star2, "path of primary star",
             verticalalignment="center", color="r")

    plt.plot([0.1*a_star2,1.1*a_star2],[1.9*a_star2,1.9*a_star2], color="g")
    plt.text(1.3*a_star2, 1.9*a_star2, "path of unseen companion",
             verticalalignment="center", color="g")

    plt.axis([-0.3*a_star2, 19.7*a_star2, 
              -3*a_star2, 3*a_star2])

    f = plt.gcf()
    f.set_size_inches(10.0,3.0)
  
    plt.savefig("astrometric_binary.png")

    
if __name__== "__main__":
    astrometric_binary()


    
        
