import numpy as np
import matplotlib.pyplot as plt

import anim_solvers.solar_system_integrator as ss

"""
Integrate the orbit of a planet around the Sun, and display the
kinetic and potential energy

In this version, we just consider orbits with the horizontal 
velocity <= the circular velocity.
"""

# we work in MKS units
G = 6.67428e-11      # m^3 kg^{-1} s^{-2}
M_sun = 1.98892e30   # kg
AU = 1.49598e11      # m
year = 3.1557e7      # s

def orbitalenergy():

    # set the semi-major axis and eccentricity
    a = 1.5874*AU
    e = 0.4

    orbit = ss.SolarSystem(GM=G*M_sun, year=year)

    orbit.add_planet(a, e, loc="perihelion")

    # compute the period of the orbit from Kepler's law and make
    # the timestep by 1/720th of a period
    P = orbit.period(0)

    print("period = ", P/year)

    sol = orbit.integrate(360, 2*P/year)

    # plotting

    iframe = 0

    for n in range(len(sol[0].t)):

        fig, ax = plt.subplots()

        # plot the Sun
        plt.scatter([0], [0], s=1600, marker=(20, 1), color="k")
        plt.scatter([0], [0], s=1500, marker=(20, 1), color="#FFFF00")

        # plot planet
        ax.plot(sol[0].x, sol[0].y, color="C0")
        ax.scatter([sol[0].x[n]], [sol[0].y[n]], s=100, color="C0")

        # compute the kinetic energy / kg
        KE = 0.5*(sol[0].vx[n]**2 + sol[0].vy[n]**2)

        # compute the potential energy / kg
        r = np.sqrt(sol[0].x[n]**2 + sol[0].y[n]**2)
        PE = - G*M_sun/r

        ax.axis([-3.5*AU, 2*AU, -3.5*AU, 2*AU])
        ax.set_aspect("equal", "datalim")

        fig.set_size_inches(7.2, 7.2)

        ax.set_title("Orbital Energy")

        ax.text(-2.5*AU, -2.2*AU,  f"KE / unit mass (J/kg): {KE:10.5e}", fontsize="13")
        ax.text(-2.5*AU, -2.5*AU, f"PE / unit mass (J/kg): {PE:10.5e}", fontsize="13")
        ax.text(-2.5*AU, -2.8*AU, f"total energy / unit mass (J/kg): {PE + KE:10.5e}", fontsize="13")

        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")

        ax.xaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))
        ax.yaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))

        fig.tight_layout()
        fig.savefig(f"orbitalenergy_{n:04d}.png")
        plt.close(fig)


if __name__ == "__main__":
    orbitalenergy()
