import math
import numpy as np
import matplotlib.pyplot as plt

import anim_solvers.solar_system_integrator as ss

# Integrate the orbit of a planet around the Sun, and display the
# kinetic and potential energy
#
# In this version, we just consider orbits with the horizontal 
# velocity <= the circular velocity.


# M. Zingale (2008-09-14)

# we work in MKS units
G = 6.67428e-11      # m^3 kg^{-1} s^{-2}
M_sun = 1.98892e30   # kg
AU = 1.49598e11      # m
year = 3.1557e7      # s

def orbitalenergy():

    # set the semi-major axis and eccentricity
    a = 1.5874*AU
    e = 0.4

    orbit = ss.SolarSystem(GM=G*M_sun, year=year)

    orbit.add_planet(a, e, loc="perihelion")

    # compute the period of the orbit from Kepler's law and make
    # the timestep by 1/720th of a period
    P = orbit.period(0)

    print "period = ", P/year

    sol = orbit.integrate(360, P/year)


    # ================================================================
    # plotting
    # ================================================================

    # plot the orbit
    iframe = 0

    # v1
    for n in range(len(sol[0].t)):

        plt.clf()

        # plot the foci
        plt.scatter([0], [0], s=250, marker=(5,1), color="k")
        plt.scatter([0], [0], s=200, marker=(5,1), color="y")

        # plot planet 
        plt.plot(sol[0].x, sol[0].y, color="r")
        plt.scatter([sol[0].x[n]], [sol[0].y[n]], s=100, color="r")


        # compute the kinetic energy / kg
        KE = 0.5*(sol[0].vx[n]**2 + sol[0].vy[n]**2)

        # compute the potential energy / kg
        r = math.sqrt(sol[0].x[n]**2 + sol[0].y[n]**2)
        PE = - G*M_sun/r
        
        plt.axis([-4*AU,2*AU,-4*AU,2*AU])
        
        f = plt.gcf()
        f.set_size_inches(7.2,7.2)

        ax = plt.gca()
        ax.set_aspect("equal", "datalim")

        plt.title("Orbital Energy")

        plt.text(-3.5*AU,-3*AU,  "KE / unit mass (J/kg): %10.5e" % (KE))
        plt.text(-3.5*AU,-3.3*AU,"PE / unit mass (J/kg): %10.5e" % (PE))
        plt.text(-3.5*AU,-3.6*AU,"total energy / unit mass (J/kg): %10.5e" % (PE + KE))

        plt.xlabel("x [m]")
        plt.ylabel("y [m]")

        ax.xaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))
        ax.yaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))

        plt.savefig("orbitalenergy_%04d.png" % n)

    
if __name__== "__main__":
    orbitalenergy()


    
        
