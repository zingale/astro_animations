# illustrate the orbital resonance in between Juputer and the asteroid belt

import math
import random
import matplotlib.pyplot as plt

import anim_solvers.solar_system_integrator as ssi


def asteroids():

    ss = ssi.SolarSystem()

    # Jupiter
    P_jupiter = 12 # years
    e = 0.0
    ss.add_planet_by_period(P_jupiter, e, loc="perihelion")

    # an asteroid, with some eccentricity but 1/2 the period
    e = 0.3
    ss.add_planet_by_period(0.5*P_jupiter, e, loc="perihelion")

    # add some background asteroids with random orbits
    a_min = 2.0
    a_max = 3.5
    n_asteroids = 50

    for n in range(n_asteroids):
        a = random.uniform(a_min, a_max)
        e = random.uniform(0.0, 0.4)
        if n % 2 == 0:
            loc="perihelion"
        else:
            loc="aphelion"
            
        ss.add_planet(a, e, loc=loc, rot="random")

    # integrate
    nsteps_per_year = 45.0
    sol = ss.integrate(nsteps_per_year, 2*P_jupiter)
    
    
    # plots
    iframe = 0
    
    for n in range(len(sol[0].t)):

        plt.clf()
        
        # plot the Sun
        plt.scatter([0], [0], s=1600, marker=(20,1), color="k")
        plt.scatter([0], [0], s=1500, marker=(20,1), color="#FFFF00")

        # plot Jupiter
        plt.plot(sol[0].x, sol[0].y, color="k")
        plt.scatter([sol[0].x[n]], [sol[0].y[n]], s=500, color="k")

        # plot our asteroid
        plt.plot(sol[1].x, sol[1].y, color="r")
        plt.scatter([sol[1].x[n]], [sol[1].y[n]], s=50, color="r")


        # and the background asteroids
        for k in range(2, 2+n_asteroids):
            plt.plot(sol[k].x, sol[k].y, color="0.85", zorder=-100)
            plt.scatter([sol[k].x[n]], [sol[k].y[n]], s=25, color="0.5", zorder=-100)


        f = plt.gcf()
        f.set_size_inches(7.2, 7.2)

        plt.axis("off")
        
        ax = plt.gca()
        ax.set_aspect("equal", "datalim")

        # if jupiter is at perihelion, the pause and annotate
        perihelion = False
        if n == 0:
            perihelion = True
        elif n > 0 and n < len(sol[0].t)-1:
            if sol[0].y[n]*sol[0].y[n+1] < 0.0 and sol[0].x[n] > 0.0:
                perihelion = True

        if perihelion:
            plt.text(0.5, 0.96,
                     "Jupiter and our asteroid are at their closest point",
                     horizontalalignment="center", transform=f.transFigure,
                     fontsize="large")

            plt.text(0.5, 0.92,
                     "The gravitational force on the asteroid is always",
                     horizontalalignment="center", transform=f.transFigure)            

            plt.text(0.5, 0.89,
                     "strongest and in the same direction here",            
                     horizontalalignment="center", transform=f.transFigure)            
            

            plt.arrow(sol[1].x[n], sol[1].y[n], 1.0, 0.0, color="r",
                      length_includes_head=True,
                      head_width = 0.2, width=0.05, overhang=-0.1)

            for k in range(150):
                plt.savefig("asteroids_{:04d}.png".format(iframe))
                iframe += 1
                 
        plt.savefig("asteroids_{:04d}.png".format(iframe))
        iframe += 1

if __name__ == "__main__":
    asteroids()
    
