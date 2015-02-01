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
    nsteps_per_year = 90.0
    sol = ss.integrate(nsteps_per_year, 2*P_jupiter)
    
    
    # plots
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
            plt.plot(sol[k].x, sol[k].y, color="0.75", zorder=-100)
            plt.scatter([sol[k].x[n]], [sol[k].y[n]], s=25, color="0.75", zorder=-100)


        f = plt.gcf()
        f.set_size_inches(7.2, 7.2)

        plt.axis("off")
        
        ax = plt.gca()
        ax.set_aspect("equal", "datalim")
                                                          
        plt.savefig("asteroids_{:04d}.png".format(n))


if __name__ == "__main__":
    asteroids()
    
