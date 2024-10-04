import numpy as np
import matplotlib.pyplot as plt
import anim_solvers.earth_orbit as earth_orbit

"""
Shoot a projective horizontally some distance above Earth at various
speeds and watch the resulting orbit

In this version, we just consider orbits with the horizontal
velocity <= the circular velocity.
"""

# M. Zingale (2008-09-14)

# we work in CGS units
G = 6.67e-8
M_E = 5.9742e27
R_E = 6.378e8


def orbit():

    # set the height from which we launch things
    hinit = 1.5 * R_E
    max_rad = 25 * R_E

    # the circular orbit will be our baseline
    orbit_circ = earth_orbit.Trajectory(GM=G*M_E, R_crash=R_E)

    vinit_circ = np.sqrt(G * M_E / hinit)
    tmax_circ = 2 * np.pi * hinit / vinit_circ

    # we want dt to be constant for all of our orbits so we can compare
    # the velocities in the animation
    dt = tmax_circ / 720.0

    orbit_circ.integrate(vinit_circ, hinit, dt, max_rad)

    # orbit 1 (0.25 v_c)
    orbit_1 = earth_orbit.Trajectory(GM=G*M_E, R_crash=R_E)
    vinit_1 = 0.25*vinit_circ

    orbit_1.integrate(vinit_1, hinit, dt, max_rad)

    # orbit 2 (0.5 v_c)
    orbit_2 = earth_orbit.Trajectory(GM=G*M_E, R_crash=R_E)
    vinit_2 = 0.5*vinit_circ

    orbit_2.integrate(vinit_2, hinit, dt, max_rad)

    # orbit 3 (0.75 v_c)
    orbit_3 = earth_orbit.Trajectory(GM=G*M_E, R_crash=R_E)
    vinit_3 = 0.75*vinit_circ

    orbit_3.integrate(vinit_3, hinit, dt, max_rad)

    # plotting

    img = plt.imread("earth.png")

    # plot the orbits one by one
    iframe = 0

    trajectories = [(orbit_1, r"$v = 0.25 v_\mathrm{circular}$", "C0"),
                    (orbit_2, r"$v = 0.5 v_\mathrm{circular}$", "C1"),
                    (orbit_3, r"$v = 0.75 v_\mathrm{circular}$", "C2"),
                    (orbit_circ, r"$v = v_\mathrm{circular}$", "C3")]

    for itraj in range(len(trajectories)):

        orbit, label, color = trajectories[itraj]

        for n in range(orbit.npts):

            fig, ax = plt.subplots()

            # plot previous trajectories (if they exist)

            for iold in range(itraj):
                orbit_old, _, color_old = trajectories[iold]
                ax.plot(orbit_old.x[0:orbit_old.npts-1],
                        orbit_old.y[0:orbit_old.npts-1],
                        color=color_old, alpha=0.33)

            # now the current

            ax.set_title(label, color=color, fontsize=20)

            # draw Earth -- we will use units that are in terms of Earth radii
            ax.imshow(img, extent=[-R_E, R_E, -R_E, R_E])

            # plot the current orbit
            ax.plot(orbit.x[0:n], orbit.y[0:n], color=color)

            ax.axis([-3*R_E, 3*R_E, -3*R_E, 2*R_E])

            ax.axis("off")
            ax.set_aspect("equal", "datalim")

            fig.set_size_inches(9.6, 7.2)

            ax.text(0.05, 0.06, "Earth image credit:",
                    transform=fig.transFigure,
                    fontsize=7, color="0.50")
            ax.text(0.05, 0.04, "NASA/Apollo 17", transform=fig.transFigure,
                    fontsize=7, color="0.50")

            # print the time
            ax.text(0.7, 0.1, f"time = {orbit.t[n]/3600:6.4f} hrs.",
                    transform=fig.transFigure)

            fig.savefig(f"achieveorbit_{iframe:04d}.png")
            plt.close(fig)

            iframe += 1


if __name__ == "__main__":
    orbit()
