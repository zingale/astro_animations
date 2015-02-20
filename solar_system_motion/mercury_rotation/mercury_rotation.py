import math
import numpy as np
import matplotlib.pyplot as plt
import anim_solvers.solar_system_integrator as ssi

# Integrate the orbit of Mercury around the Sun and show its rotation
# by drawing a point on the surface.


def orbitalenergy():

    # set the semi-major axis and eccentricity
    a = 0.387
    e = 0.2056

    ss = ssi.SolarSystem()
    ss.add_planet(a, e, loc="perihelion")

    P_orbital = ss.period(0)

    # set the rotation period
    P_rotation = (2./3.)*P_orbital

    omega = 2*math.pi/P_rotation

    print "orbital period = ", P_orbital
    
    nsteps_per_year = int(360/P_orbital)
    nyears = 4*P_orbital
    sol = ss.integrate(nsteps_per_year, nyears)


    # ================================================================
    # plotting
    # ================================================================

    for n in range(len(sol[0].t)):

        plt.clf()

        # plot the foci
        plt.scatter([0], [0], s=1600, marker=(20,1), color="k")
        plt.scatter([0], [0], s=1500, marker=(20,1), color="#FFFF00")

        # plot the orbit
        plt.plot(sol[0].x, sol[0].y, color="0.5")

        # plot planet
        theta = np.radians(np.arange(360))

        r = 0.03  # exaggerate the planet's size
        x_surface = sol[0].x[n] + r*np.cos(theta)
        y_surface = sol[0].y[n] + r*np.sin(theta)
        plt.fill(x_surface, y_surface, color="r", alpha=1.0, zorder=1000)

        # plot a point on the planet's surface
        xpt = sol[0].x[n] + r*np.cos(omega*sol[0].t[n]+math.pi)
        ypt = sol[0].y[n] + r*np.sin(omega*sol[0].t[n]+math.pi)
        plt.scatter([xpt], [ypt],s=25, color="k")

        plt.axis([-0.5, 0.5, -0.5, 0.5])

        plt.axis("off")

        ax = plt.gca()
        ax.set_aspect("equal", "datalim")

        plt.subplots_adjust(left=0.05, right=0.98, bottom=0.05, top=0.98)
        
        f = plt.gcf()

        plt.text(0.5, 0.925, r"Mercury: $P_\mathrm{rotation} = (2/3) P_\mathrm{orbital}$",
                 transform=f.transFigure, horizontalalignment="center")

        f.set_size_inches(7.2,7.2)

        plt.savefig("mercury_rotation_{:04d}.png".format(n))


if __name__== "__main__":
    orbitalenergy()
