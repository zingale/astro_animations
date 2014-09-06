import math
import numpy
import matplotlib.pyplot as plt

import anim_solvers.binary_integrator as bi

# compute the orbits of a planet and the parent star. 
#
# Here we put the center of mass at the origin.  We exaggerate the
# the mass of the planet to make the star's motion more pronounced.

# M. Zingale (2008-11-30)

# we work in CGS units
G = bi.G
M_sun = bi.M_sun
AU = bi.AU
year = bi.year


def radial_velocity():

    # set the masses
    M_star = M_sun           # star's mass
    M_p    = 0.15*M_sun       # planet's mass

    # set the semi-major axis of the planet (and derive that of the star)
    a_p = 1.0*AU
    a_star = (M_p/M_star)*a_p  # M_p a_p = -M_star a_star (center of mass)

    ecc = 0.0

    theta = 0.0

    # create the solar system container
    b = bi.Binary(M_star, M_p, a_star + a_p, ecc, theta)

    
    # set the timestep in terms of the orbital period
    dt = b.P/360.0        
    tmax = 2.0*b.P  # maximum integration time

    s, p = b.integrate(dt, tmax)


    # ================================================================
    # plotting
    # ================================================================

    iframe = 0

    # first plot a legend intro sequence
    while iframe < 100:

        plt.clf()
        
        plt.subplot(121)

        plt.subplots_adjust(left=0.1,right=0.95,bottom=0.15,top=0.85,
                            wspace=0.45,hspace=0.0)

        a = plt.gca()
        a.set_aspect("equal", "datalim")
        plt.axis("off")

        plt.scatter([0], [0], s=150, marker="x", color="k")


        # plot the planet's orbit and position
        plt.plot(p.x, p.y, color="g")
        plt.scatter([p.x[0]], [p.y[0]], s=100, color="g")

        # plot the star's orbit and position
        plt.plot(s.x, s.y, color="r")
        plt.scatter([s.x[0]], [s.y[0]], s=200, color="r")

        # plot a reference line
        plt.plot([0,1*AU], [-1.2*a_p, -1.2*a_p], color="k")
        plt.text(0.5*AU, -1.4*a_p, "1 AU", horizontalalignment='center')

        plt.arrow(0.0, 1.5*a_p, 0.0, 0.2*a_p, width=0.025*a_p, 
                  facecolor="b", edgecolor="k", alpha=1.0, clip_on=0,
                  length_includes_head=True, head_width=0.075*a_p, head_length=0.075*a_p)
        plt.text(0.0, 1.4*a_p, "observer", horizontalalignment='center')

        plt.text(-1.8*a_p, -1.6*a_p, "time = %6.3f yr" % (p.t[0]/year))

        plt.axis([-1.5*a_p, 1.5*a_p, -1.5*a_p, 1.5*a_p])


        plt.subplot(122)
        
        plt.scatter([0.1], [0.8], s=200, color="r")
        plt.text(0.2, 0.8, "star", color="r")

        plt.scatter([0.1], [0.6], s=100, color="g")
        plt.text(0.2, 0.6, "planet", color="g")

        plt.text(0.2, 0.5, "(note: planet's mass exaggerated", size="small", color="k")
        plt.text(0.2, 0.45, " to enhance effect)", size="small", color="k")


        plt.scatter([0.1], [0.3], s=150, marker="x", color="k")
        plt.text(0.2, 0.3, "center of mass", color="k")

        plt.axis([0,1, 0,1])
        plt.axis("off")

        f = plt.gcf()
        f.set_size_inches(12.8,7.2)

        plt.savefig("radial_velocity_%04d.png" % iframe)

        iframe += 1


    # now plot the actual data
    for n in range(len(s.t)):

        plt.clf()

        plt.subplot(121)

        plt.subplots_adjust(left=0.1,right=0.95,bottom=0.15,top=0.85,
                            wspace=0.45,hspace=0.0)

        a = plt.gca()
        a.set_aspect("equal", "datalim")
        plt.axis("off")

        plt.scatter([0], [0], s=150, marker="x", color="k")


        # plot the planet's orbit and position
        plt.plot(p.x, p.y, color="g")
        plt.scatter([p.x[n]], [p.y[n]], s=100, color="g")

        # plot the star's orbit and position
        plt.plot(s.x, s.y, color="r")
        plt.scatter([s.x[n]], [s.y[n]], s=200, color="r")

        # plot a reference line
        plt.plot([0,1*AU], [-1.2*a_p,-1.2*a_p], color="k")
        plt.text(0.5*AU, -1.4*a_p, "1 AU", horizontalalignment='center')

        plt.arrow(0.0, 1.5*a_p, 0.0, 0.2*a_p, width=0.025*a_p, 
                  facecolor="b", edgecolor="k", alpha=1.0, clip_on=0,
                  length_includes_head=True, head_width=0.075*a_p, head_length=0.075*a_p)
        plt.text(0.0, 1.4*a_p, "observer", horizontalalignment='center')

        plt.text(-1.8*a_p, -1.6*a_p, "time = %6.3f yr" % (s.t[n]/year))

        plt.axis([-1.5*a_p, 1.5*a_p, -1.5*a_p, 1.5*a_p])


        plt.subplot(122)

        # plot km/s vs years
        # note: radial velocity convention is that if it is moving 
        # toward us, then the velocity is negative.  Since the observer
        # is at +Y, we want to plot -vy as the radial velocity.
        plt.plot(s.t/b.P,-s.vy/1.e5, color="r")
        plt.scatter([s.t[n]/b.P], [-s.vy[n]/1.e5], color="r", s=50)

        # plot a reference line at vy = 0
        plt.plot([0.0,tmax/b.P], [0.0,0.0], "k--")


        plt.axis([0.0, tmax/b.P,
                    int(numpy.min(s.vy))/1.e5 - 1,
                    int(numpy.max(s.vy))/1.e5 + 1])


        plt.xlabel("t/period")
        plt.ylabel("radial velocity [km/s]")

        f = plt.gcf()
        f.set_size_inches(12.8,7.2)
        
        plt.savefig("radial_velocity_%04d.png" % iframe)

        iframe += 1

    
if __name__== "__main__":
    radial_velocity()


    
        
