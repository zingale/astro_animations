import numpy as np
import matplotlib.pyplot as plt
import anim_solvers.earth_orbit as earth_orbit

"""
Shoot a projective horizontally some distance above Earth at various
speeds and watch the resulting orbit

Give it a boost a perihelion and observe the new orbit.
"""

# we work in CGS units
G = 6.67e-8
M_E = 5.9742e27
R_E = 6.378e8
SMALL = 1.e-12

max_coord = 3*R_E

def orbit():

    # set the height from which we launch things
    hinit = 1.5 * R_E
    max_rad = 25 * R_E

    # the circular orbit will be our baseline
    vinit_circ = np.sqrt(G * M_E / hinit)
    tmax_circ = 2 * np.pi * hinit / vinit_circ

    # we want dt to be constant for all of our orbits so we can compare
    # the velocities in the animation
    dt = tmax_circ / 240.0

    # orbit 1: v = v_c
    orbit1 = earth_orbit.Trajectory(GM=G*M_E, R_crash=R_E)
    vinit1 = vinit_circ

    orbit1.integrate(vinit1, hinit, dt, max_rad)

    # orbit 2: v = 1.1*v_c
    orbit2 = earth_orbit.Trajectory(GM=G*M_E, R_crash=R_E)
    vinit2 = 1.1*vinit_circ

    orbit2.integrate(vinit2, hinit, dt, max_rad)

    # orbit 3: v = (1.1)**2*v_c
    orbit3 = earth_orbit.Trajectory(GM=G*M_E, R_crash=R_E)
    vinit3 = 1.1*1.1*vinit_circ

    orbit3.integrate(vinit3, hinit, dt, max_rad)

    # plotting

    img = plt.imread("earth.png")

    # plot the orbits one by one
    iframe = 0

    trajectories = [(orbit1, r"$\mathrm{circular\ orbit:}\ v_\mathrm{peri} = v_\mathrm{circular}$", "C0"),
                    (orbit2, r"$\mathrm{first\ boost:}\ v_\mathrm{peri} = 1.1 \times\ v_\mathrm{circular}$", "C1"),
                    (orbit3, r"$\mathrm{second\ boost:}\ v_\mathrm{peri} = (1.1)^2 \times\ v_\mathrm{circular}$", "C2")]

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
            ax.plot(orbit.x[0:n+1], orbit.y[0:n+1], color=color)

            # draw the spaceship
            ax.scatter([orbit.x[n]], [orbit.y[n]], color="b", marker="o")

            angle = np.degrees(np.arctan2(orbit.y[n], (orbit.x[n] + SMALL)))

            # draw flames for the boost at start of this orbit
            if itraj > 0:
                if (78.5 < angle <= 90 and n < 0.25*orbit2.npts):
                    ax.scatter([orbit.x[n]-0.025*max_coord], [orbit.y[n]],
                               color="r", marker="<")

            if itraj < len(trajectories)-1:
                # draw flames for the boost at the start of next orbit
                if (89 <= angle < 102.5 and n > 0.75*orbit.npts):
                    ax.scatter([orbit.x[n]-0.025*max_coord], [orbit.y[n]],
                               color="r", marker="<")

            ax.axis("off")
            ax.set_aspect("equal")

            fig.set_size_inches(9.6, 7.2)

            ax.set_xlim(-3*R_E, 3*R_E)
            ax.set_ylim(-4.5*R_E, 2.5*R_E)

            # print the image credit
            ax.text(0.05, 0.06, "Earth image credit:",
                    transform=fig.transFigure,
                    fontsize=7, color="0.50")
            ax.text(0.05, 0.04, "NASA/Apollo 17", transform=fig.transFigure,
                    fontsize=7, color="0.50")

            # print the time
            ax.text(0.7, 0.1, f"time = {orbit.t[n]/3600:6.4f} hrs.",
                    transform=fig.transFigure)

            fig.savefig(f"changeorbit_{iframe:04d}.png")
            plt.close(fig)

            iframe += 1


if __name__ == "__main__":
    orbit()
