#!/bin/env python

"""compute the orbit of a 2 planets around the Sun using, showing the
difference in properties with semi-major axis and eccentricity

"""


import matplotlib.pyplot as plt

import anim_solvers.solar_system_integrator as solar_system_integrator


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

    xmin = min(sol[0].x.min(), sol[1].x.min())
    xmax = max(sol[0].x.max(), sol[1].x.max())

    ymin = min(sol[0].y.min(), sol[1].y.min())
    ymax = max(sol[0].y.max(), sol[1].y.max())

    # plotting
    for n in range(len(sol[0].x)):

        fig, ax = plt.subplots()

        # plot the foci
        ax.scatter([0], [0], s=250, marker=(5, 1), color="k")
        ax.scatter([0], [0], s=200, marker=(5, 1), color="y")

        # plot planet A
        ax.plot(sol[0].x, sol[0].y, color="C0")
        ax.scatter([sol[0].x[n]], [sol[0].y[n]], s=100, color="C0")

        # plot planet B
        ax.plot(sol[1].x, sol[1].y, color="C1")
        ax.scatter([sol[1].x[n]], [sol[1].y[n]], s=100, color="C1")

        ax.axis([1.1 * xmin, 1.1 * xmax, 1.1 * ymin, 1.1 * ymax])
        ax.axis("off")
        ax.set_aspect("equal", "datalim")

        fig.set_size_inches(9.6, 7.2)

        ax.text(0.05, 0.05, f"time = {sol[0].t[n]:6.3f} yr",
                transform=fig.transFigure)
        ax.text(0.05, 0.9, f"a = {a_A:6.3f}, e = {ecc_A:5.2f}", color="C0",
                transform=fig.transFigure)
        ax.text(0.05, 0.86, f"a = {a_B:6.3f}, e = {ecc_B:5.2f}", color="C1",
                transform=fig.transFigure)

        fig.savefig(f"orbit_{n:04d}.png")
        plt.subplots_adjust(left=0.025, right=0.975, bottom=0.025, top=0.975)
        plt.close(fig)


if __name__ == "__main__":
    doit()
