# illustrate the orbital resonance in between Juputer and the asteroid belt

import random
import matplotlib.pyplot as plt

import anim_solvers.solar_system_integrator as ssi


def asteroids():

    ss = ssi.SolarSystem()

    # Jupiter
    P_jupiter = 12  # years
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
            loc = "perihelion"
        else:
            loc = "aphelion"

        ss.add_planet(a, e, loc=loc, rot="random")

    # integrate
    nsteps_per_year = 45.0
    sol = ss.integrate(nsteps_per_year, 4*P_jupiter)

    # plots
    iframe = 0

    for n in range(len(sol[0].t)):

        fig, ax = plt.subplots()

        # plot the Sun
        plt.scatter([0], [0], s=1600, marker=(20, 1), color="k")
        plt.scatter([0], [0], s=1500, marker=(20, 1), color="#FFFF00")

        # plot Jupiter
        plt.plot(sol[0].x, sol[0].y, color="C1")
        plt.scatter([sol[0].x[n]], [sol[0].y[n]], s=500, color="C1")

        # plot our asteroid
        plt.plot(sol[1].x, sol[1].y, color="C0")
        plt.scatter([sol[1].x[n]], [sol[1].y[n]], s=50, color="C0")

        # and the background asteroids
        for k in range(2, 2+n_asteroids):
            plt.plot(sol[k].x, sol[k].y, color="0.85", zorder=-100)
            plt.scatter([sol[k].x[n]], [sol[k].y[n]], s=25,
                        color="0.5")

        fig.set_size_inches(7.2, 7.2)

        ax.axis("off")
        ax.set_aspect("equal", "datalim")

        # if jupiter is at perihelion, the pause and annotate
        perihelion = False
        if n == 0:
            perihelion = True
        elif n > 0 and n < len(sol[0].t)-1:
            if sol[0].y[n]*sol[0].y[n+1] < 0.0 and sol[0].x[n] > 0.0:
                perihelion = True

        fig.subplots_adjust(left=0.0125, right=0.985, bottom=0.0125, top=0.985)

        if perihelion:
            ax.text(0.5, 0.97,
                    "Jupiter and our asteroid are at their closest point",
                    horizontalalignment="center", transform=fig.transFigure,
                    fontsize="large", zorder=1000)

            ax.text(0.5, 0.94,
                    "The gravitational force on the asteroid is strongest",
                    horizontalalignment="center",
                    transform=fig.transFigure, zorder=1000)

            ax.arrow(sol[1].x[n], sol[1].y[n], 1.0, 0.0, color="C0",
                     length_includes_head=True,
                     head_width=0.2, width=0.05, overhang=-0.1)

            for k in range(150):
                fig.savefig(f"asteroids_{iframe:04d}.png")
                iframe += 1

        fig.savefig(f"asteroids_{iframe:04d}.png")
        plt.close(fig)
        iframe += 1


if __name__ == "__main__":
    asteroids()
