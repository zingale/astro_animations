import math
import numpy as np
import matplotlib.pyplot as plt

import anim_solvers.solar_system_integrator as ssi

# Show the orbit of earth with its axis around the Sun

def seasons():

    # set the semi-major axis and eccentricity
    a = 1.0  # AU
    e = 0.0

    s = ssi.SolarSystem()

    s.add_planet(a, e, loc="perihelion")

    P = s.period(0)

    nsteps_per_year = 360
    num_years = 1
    
    sol = s.integrate(nsteps_per_year, num_years)
    
    # apply a projection to account for the inclination
    # wrt the observer
    inc = 80
    
    sol[0].y = sol[0].y*np.cos(np.radians(inc))
        
    # ================================================================
    # plotting
    # ================================================================

    # plot the orbit
    for n in range(len(sol[0].t)):

        plt.clf()

        # plot the Sun at the foci
        plt.scatter([0], [0], s=2600, marker=(20,1), color="k", zorder=0)
        plt.scatter([0], [0], s=2500, marker=(20,1), color="#FFFF00", zorder=0)

        # plot the orbit
        plt.plot(sol[0].x, sol[0].y, color="0.5", ls="--", zorder=-50)


        # plot planet -- hide it with zorder
        if sol[0].y[n] > 0:
            z = -20
        else:
            z = 20
            
        theta = np.radians(np.arange(360))
        r = 0.075  # exaggerate the planet's size
        x_surface = sol[0].x[n] + r*np.cos(theta)
        y_surface = sol[0].y[n] + r*np.sin(theta)
        plt.fill(x_surface, y_surface,"b", edgecolor="b", zorder=z)

        # axis
        tilt = np.radians(23.5)
        L = 0.1
        x = [-L*np.sin(tilt), L*np.sin(tilt)]
        y = [-L*np.cos(tilt), L*np.cos(tilt)]

        plt.plot(sol[0].x[n]+x, sol[0].y[n]+y, color="k", lw=2, zorder=z)

        # equator
        y = [ r*np.sin(tilt), -r*np.sin(tilt)]
        x = [-r*np.cos(tilt),  r*np.cos(tilt)]

        plt.plot(sol[0].x[n]+x, sol[0].y[n]+y, color="k", lw=1, ls="--", zorder=z)        
        
        plt.axis([-1.2, 1.2, -0.9, 0.9])

        f = plt.gcf()
        f.set_size_inches(9.6, 7.2)

        plt.title("Seasons")

        plt.axis("off")

        ax = plt.gca()
        ax.set_aspect("equal", "datalim")

        plt.savefig("seasons_{:04d}".format(n))    

    
if __name__== "__main__":
    seasons()


    
        
