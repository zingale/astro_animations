#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import anim_solvers.solar_system_integrator as ssi

# demonstrate synchronous rotation of the Moon

# we work in MKS units
d = 3.84399e10           # semi-major axis [cm]
AU = 1.5e13              # cm
M_earth = 5.972e27       # g
M_sun = 2.e33

def doit():

    # set the semi-major axis and eccentricity
    a = d / AU
    e = 0.0

    # SolarSystem assumes orbits around the Sun so we need to rescale
    # it so we work in Earth units

    ss = ssi.SolarSystem(GM=4*np.pi**2*M_earth/M_sun)
    ss.add_planet(a, e, loc="perihelion")

    # compute the period of the orbit from Kepler's law and make
    # the timestep by 1/720th of a period
    P_orbital = ss.period(0)

    # set the rotation period
    P_rotation = P_orbital

    omega = 2*np.pi/P_rotation

    print("period = ", P_orbital * 365)

    nsteps_per_period = 720
    tmax = 2 * P_orbital

    nsteps_per_year = int(nsteps_per_period / P_orbital)

    sol = ss.integrate(nsteps_per_year, tmax)

    # plotting

    for n in range(len(sol[0].t)):

        fig, ax = plt.subplots()

        # plot the Earth
        ax.scatter([0], [0], s=650, color="C0")

        # plot the orbit
        ax.plot(sol[0].x * AU, sol[0].y * AU,
                color="0.5", ls=":", lw=1)

        # plot moon (use zorder to put this on top of the orbit line)
        theta = np.arange(180)
        r = 0.05*d  # exaggerate the moon's size
        x_surface = sol[0].x[n] * AU + r*np.cos(theta)
        y_surface = sol[0].y[n] * AU + r*np.sin(theta)
        ax.fill(x_surface, y_surface, "0.7",
                edgecolor="0.7", alpha=1.0, zorder=1000)

        # plot a point on the moon's surface
        xpt = sol[0].x[n] * AU + r*np.cos(omega*sol[0].t[n] + np.pi)
        ypt = sol[0].y[n] * AU + r*np.sin(omega*sol[0].t[n] + np.pi)
        ax.scatter([xpt], [ypt], s=25, color="k", zorder=10000)

        ax.axis([-1.1*d, 1.1*d, -1.1*d, 1.1*d])
        ax.axis("off")

        fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)

        ax.set_aspect("equal", "datalim")

        fig.set_size_inches(7.2, 7.2)

        fig.savefig(f"moon_rotation_{n:04d}.png")
        plt.close(fig)


if __name__ == "__main__":
    doit()
