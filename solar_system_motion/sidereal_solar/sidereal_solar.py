import math
import numpy as np
import matplotlib.pylab as plt
import anim_solvers.solar_system_integrator as ssi

# consider a planet whose sidereal rotation period = 1/2 its orbital
# period.  Its solar day would then be 1 orbital period.  This
# animation demonstrates that.


def sidereal_solar():

    # set the semi-major axis and eccentricity
    a = 1.0
    e = 0.0

    ss = ssi.SolarSystem()
    ss.add_planet(a, e, loc="perihelion")

    
    # compute the period of the orbit from Kepler's law and ...
    P_orbital = ss.period(0)

    # ... set the rotation period
    P_rotation = 0.5*P_orbital

    omega = 2*math.pi/P_rotation


    # integrate
    nsteps_per_year = 720
    num_years = 2
    
    sol = ss.integrate(nsteps_per_year, num_years)


    # ================================================================
    # plotting
    # ================================================================


    for n in range(len(sol[0].t)):

        plt.clf()

        # plot the foci
        plt.scatter([0], [0], s=250, marker=(5,1), color="k")
        plt.scatter([0], [0], s=200, marker=(5,1), color="y")

        # plot the orbit
        plt.plot(sol[0].x, sol[0].y, color="0.5", ls="--")

        # plot planet 
        theta = np.radians(np.arange(360))                
        r = 0.05  # exaggerate the planet's size
        x_surface = sol[0].x[n] + r*np.cos(theta)
        y_surface = sol[0].y[n] + r*np.sin(theta)
        plt.fill(x_surface, y_surface, "c", edgecolor="c", zorder=100)

        # plot a point on the planet's surface
        xpt = sol[0].x[n] + r*np.cos(omega*sol[0].t[n]+math.pi/2.0)
        ypt = sol[0].y[n] + r*np.sin(omega*sol[0].t[n]+math.pi/2.0)

        plt.scatter([xpt], [ypt], s=25, color="k", zorder=100)

        plt.axis([-1.2, 1.2, -1.2, 1.2])

        ax = plt.gca()
        ax.set_aspect("equal", "datalim")
        
        plt.subplots_adjust(left=0.05,right=0.98,bottom=0.05,top=0.98)     
        
        f = plt.gcf()
        f.set_size_inches(7.2, 7.2)

        plt.axis("off")

        plt.savefig("sidereal_solar_{:04d}.png".format(n))

    
if __name__== "__main__":
    sidereal_solar()


    
        
