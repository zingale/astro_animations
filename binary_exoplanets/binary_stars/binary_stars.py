import math
import numpy as np
import pylab

import anim_solvers.binary_integrator as bi

M_sun = bi.M_sun
AU = bi.AU
year = bi.year
G = bi.G

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

def find_scinotat(number):

    b = int(math.log10(math.fabs(number)))
    a = number/10**b
    
    return a, b


def radial_velocity():

    # set the masses
    M_star1 = 2.0*M_sun      # star 1's mass
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
    annotate = False

    # create the binary object
    b = bi.Binary(M_star1, M_star2, a_star1 + a_star2, ecc, theta, annotate=annotate)


    # set the timestep in terms of the orbital period
    dt = b.P/360.0        
    tmax = 2.0*b.P  # maximum integration time

    s1, s2 = b.integrate(dt, tmax)


    # ================================================================
    # plotting
    # ================================================================

    pylab.rc('text', usetex=True)
    pylab.rc('font', family='serif')

    iframe = 0

    for n in range(len(s1.t)):

        pylab.clf()

        pylab.subplots_adjust(left=0.1,right=0.9,bottom=0.1,top=0.9)

        a = pylab.gca()
        a.set_aspect("equal", "datalim")
        pylab.axis("off")

        pylab.scatter([0], [0], s=150, marker="x", color="k")

        # if e = 0 and M_star1 = M_star2, then the orbits lie on top of one
        # another, so plot only a single orbital line.

        # plot star 1's orbit and position
        symsize = 200
        if not (b.M1 == b.M2 and b.e == 0.0):
            pylab.plot(s1.x, s1.y, color="r")
        else:
            pylab.plot(s1.x, s1.y, color="k")

        pylab.scatter([s1.x[n]], [s1.y[n]], s=symsize, color="r")

        # plot star 2's orbit and position
        symsize = 200*(b.M2/b.M1)
        if not (b.M1 == b.M2 and b.e == 0.0):
            pylab.plot(s2.x, s2.y, color="g")

        pylab.scatter([s2.x[n]], [s2.y[n]], s=symsize, color="g")


        if annotate:
            # plot a reference line
            pylab.plot([0,1*AU], [-1.2*b.a2,-1.2*b.a2], color="k")
            pylab.text(0.5*AU, -1.4*b.a2, "1 AU", 
                       horizontalalignment='center')

            # display time
            pylab.text(-1.4*b.a2, -1.3*b.a2, 
                        "time = %6.3f yr" % (s1.t[n]/year))

        # display information about stars
        pylab.text(-1.4*b.a2, 1.3*b.a2, 
                    r"mass ratio: %3.2f" % (b.M1/b.M2), color="k")
        pylab.text(-1.4*b.a2, 1.1*b.a2, 
                    r"eccentricity: %3.2f" % (b.e), color="k")


        # energies
        if annotate:
            KE1 = 0.5*b.M1*(s1.vx[n]**2 + s1.vy[n]**2)
            KE2 = 0.5*b.M2*(s2.vx[n]**2 + s2.vy[n]**2)
            PE = -G*b.M1*b.M2/ \
                math.sqrt((s1.x[n] - s2.x[n])**2 + (s1.y[n] - s2.y[n])**2)

            print KE1, KE2, PE, KE1 + KE2 + PE

            # KE 1
            sig, ex = find_scinotat(KE1)
            pylab.text(0, 1.3*b.a2, r"$K_1 =$")
            pylab.text(0.3*b.a2, 1.3*b.a2, r"$%+4.2f \times 10^{%2d}$ erg" % (sig,ex))

            sig, ex = find_scinotat(KE2)
            pylab.text(0, 1.15*b.a2, r"$K_2 =$")
            pylab.text(0.3*b.a2, 1.15*b.a2, r"$%+4.2f \times 10^{%2d}$ erg" % (sig,ex))

            sig, ex = find_scinotat(PE)            
            pylab.text(0, 1.0*b.a2, r"$U =$")
            pylab.text(0.3*b.a2, 1.0*b.a2, r"$%+4.2f \times 10^{%2d}$ erg" % (sig,ex))

            sig, ex = find_scinotat(KE1 + KE2 + PE)            
            pylab.text(0, 0.85*b.a2, r"$E =$")
            pylab.text(0.3*b.a2, 0.85*b.a2, r"$%+4.2f \times 10^{%2d}$ erg" % (sig,ex))

        pylab.axis([-1.4*b.a2,1.4*b.a2,-1.4*b.a2,1.4*b.a2])

        f = pylab.gcf()
        f.set_size_inches(7.2,7.2)
  
        pylab.savefig("binary_star_%04d.png" % iframe)

        iframe += 1

    
if __name__== "__main__":
    radial_velocity()


    
        
