import math
import numpy as np
import matplotlib.pyplot as plt

import anim_solvers.binary_integrator as bi

# compute the orbits of stars in a binary system.
#
# Here we put the star 1 at the origin -- this shows what the orbit of
# the second star looks like from the point of view of the first.
#
# This version allows for elliptical orbits with some arbitrary
# orientation wrt to the observer (although, still face-on)
#
# The plotting assumes that M_2 < M_1
#  
# M. Zingale (2009-02-12)

# we work in CGS units

M_sun = bi.M_sun
AU = bi.AU
year = bi.year
G = bi.G


# a simple class to serve as a container for the orbital information
# for the two stars

def find_scinotat(number):

    b = int(math.log10(math.fabs(number)))
    a = number/10**b
    
    return a, b


def radial_velocity():

    # set the masses
    M_star1 = 4.0*M_sun      # star 1's mass
    M_star2 = M_sun      # star 2's mass

    # set the semi-major axis of the star 2 (and derive that of star 1)
    # M_star2 a_star2 = -M_star1 a_star1 (center of mass)
    a_star2 = 10.0*AU
    a_star1 = (M_star2/M_star1)*a_star2  

    # set the eccentricity
    ecc = 0.4

    # set the angle to rotate the semi-major axis wrt the observer
    theta = math.pi/6.0

    # display additional information
    annotate = True

    # create the binary object
    b = bi.Binary(M_star1, M_star2, a_star1 + a_star2, ecc, theta, 
                  annotate=annotate)


    
    # set the timestep in terms of the orbital period
    dt = b.P/360.0        
    tmax = 2.0*b.P  # maximum integration time

    s1, s2 = b.integrate(dt, tmax)


    # ================================================================
    # plotting
    # ================================================================

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')


    iframe = 0

    for n in range(len(s1.t)):
        
        plt.clf()

        plt.subplots_adjust(left=0.1,right=0.9,bottom=0.1,top=0.9)

        a = plt.gca()
        a.set_aspect("equal", "datalim")
        plt.axis("off")


        # offsets -- we want star 1 at the origin
        xoffset = s1.x[n]
        yoffset = s1.y[n]

        plt.scatter([0-xoffset], [0-yoffset], s=150, marker="x", color="k")

        # plot star 1's position
        symsize = 200
        plt.scatter([s1.x[n]-xoffset], [s1.y[n]-yoffset], s=symsize, color="r")

        # plot star 2's orbit and position
        symsize = 200*(b.M2/b.M1)
        plt.plot(s2.x-s1.x, s2.y-s1.y, color="g")

        plt.scatter([s2.x[n]-xoffset], [s2.y[n]-yoffset], s=symsize, color="g")

        if annotate:
            # plot a reference line
            plt.plot([0,1*AU], [-1.6*b.a2,-1.6*b.a2], color="k")
            plt.text(0.5*AU, -1.8*b.a2, "1 AU", 
                     horizontalalignment='center')

            # display time
            plt.text(-1.8*b.a2, -1.7*b.a2, 
                        "time = %6.3f yr" % (s1.t[n]/year))

        # display information about stars
        plt.text(-1.8*b.a2, 0.9*b.a2, 
                    r"mass ratio: %3.2f" % (b.M1/b.M2), color="k")
        plt.text(-1.8*b.a2, 0.7*b.a2, 
                    r"eccentricity: %3.2f" % (b.e), color="k")


        plt.axis([-1.8*b.a2, 1.0*b.a2, -1.8*b.a2, 1.0*b.a2])

        f = plt.gcf()
        f.set_size_inches(7.2,7.2)
  
        plt.savefig("binary_star_%04d.png" % iframe)

        iframe += 1

    
if __name__== "__main__":
    radial_velocity()


    
        
