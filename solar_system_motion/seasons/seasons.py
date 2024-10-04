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

    # plotting

    # plot the orbit
    for n in range(len(sol[0].t)):

        fig, ax = plt.subplots()

        # plot the Sun at the foci
        ax.scatter([0], [0], s=3500, marker=(20, 1), color="k", zorder=0)
        ax.scatter([0], [0], s=3400, marker=(20, 1), color="#FFFF00", zorder=0)

        # plot the orbit
        ax.plot(sol[0].x, sol[0].y, color="0.5", ls="--", zorder=-50)

        # plot planet -- hide it with zorder
        if sol[0].y[n] > 0:
            z = -20
        else:
            z = 20

        theta = np.radians(np.arange(360))
        r = 0.05  # exaggerate the planet's size
        x_surface = sol[0].x[n] + r*np.cos(theta)
        y_surface = sol[0].y[n] + r*np.sin(theta)
        ax.fill(x_surface, y_surface, "C0", edgecolor="C0", zorder=z)

        # axis
        tilt = np.radians(23.5)
        L = 0.1
        x = [-L*np.sin(tilt), L*np.sin(tilt)]
        y = [-L*np.cos(tilt), L*np.cos(tilt)]

        ax.plot(sol[0].x[n]+x, sol[0].y[n]+y, color="k", lw=2, zorder=z)

        # equator
        y = [ r*np.sin(tilt), -r*np.sin(tilt)]
        x = [-r*np.cos(tilt),  r*np.cos(tilt)]

        ax.plot(sol[0].x[n]+x, sol[0].y[n]+y, color="k", lw=1, ls="--", zorder=z)

        ax.axis([-1.2, 1.2, -0.9, 0.9])
        ax.axis("off")
        ax.set_aspect("equal", "datalim")

        fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)

        fig.set_size_inches(12.8, 7.2)

        ax.set_title("Seasons", fontsize="14")

        fig.savefig("seasons_{:04d}".format(n))
        plt.close(fig)


if __name__ == "__main__":
    seasons()
