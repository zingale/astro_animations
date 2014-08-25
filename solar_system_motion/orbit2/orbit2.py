#!/bin/env python

import math
import numpy
import pylab

import anim_solvers.solar_system_integrator as solar_system_integrator

# compute the orbit of a 2 planets around the Sun using, showing the difference in 
# properties with semi-major axis and eccentricity

# M. Zingale


def doit():

    # planet data
    ecc_A = 0.0   # eccentricity of planet A
    ecc_B = 0.4   # eccecntricity of planet B
    
    a_A = 1.0     # semi-major axis of planet A
    a_B = 4.0**(1./3.)  # semi-major axis of planet B

    # integration data
    nsteps_year = 365   # number of steps per year
    nyears = 4          # total integration time (years)

    s = solar_system_integrator.SolarSystem()

    s.add_planet(a_A, ecc_A, loc="perihelion")
    s.add_planet(a_B, ecc_B, loc="aphelion")

    sol = s.integrate(nsteps_year, nyears)



    # plotting
    for n in range(len(sol[0].x)):

        pylab.clf()

        # plot the foci
        pylab.scatter([0], [0], s=250, marker=(5,1), color="k")
        pylab.scatter([0], [0], s=200, marker=(5,1), color="y")

        # plot planet A
        pylab.plot(sol[0].x, sol[0].y, color="r")
        pylab.scatter([sol[0].x[n]], [sol[0].y[n]], s=100, color="r")

        # plot planet B
        pylab.plot(sol[1].x, sol[1].y, color="b")
        pylab.scatter([sol[1].x[n]], [sol[1].y[n]], s=100, color="b")

        pylab.axis([-2.5,1.5,-1.8,1.8])

        pylab.axis("off")
        ax = pylab.gca()
        ax.set_aspect("equal", "datalim")

        f = pylab.gcf()
        f.set_size_inches(9.6,7.2)

        pylab.text(0.05, 0.05, "time = %6.3f yr" % sol[0].t[n], transform=f.transFigure)
        pylab.text(0.05, 0.9, "a = %6.3f, e = %5.2f" % (a_A, ecc_A), color="r", 
                   transform=f.transFigure)
        pylab.text(0.05, 0.85, "a = %6.3f, e = %5.2f" % (a_B, ecc_B), color="b", 
                   transform=f.transFigure)

        pylab.savefig("orbit_%04d.png" % n)

    
if __name__== "__main__":
    doit()


    
        
