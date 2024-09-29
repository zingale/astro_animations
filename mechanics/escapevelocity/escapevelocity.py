import numpy as np
import matplotlib.pyplot as plt
import anim_solvers.earth_orbit as earth_orbit

"""
Shoot a projective horizontally some distance above Earth at various
speeds and watch the resulting orbit

Important speeds are the circular velocity (v_circ) and the escape
velocity (v_escp)

We consider:
             v < v_circ   (crashes into Earth)
             v = v_circ   (perfectly circular orbit)
    v_circ < v < v_escp   (elliptical orbit)
             v > v_escp   (escapes Earth's gravity)
"""

# we work in CGS units
G = 6.67e-8
M_E = 5.9742e27
R_E = 6.378e8


def orbit_plot():

    # set the height from which we launch things
    hinit = 1.5 * R_E
    max_rad = 25 * R_E

    # the circular orbit will be our baseline
    orbit_circ = earth_orbit.Trajectory(GM=G*M_E, R_crash=R_E)

    vinit_circ = np.sqrt(G * M_E / hinit)
    tmax_circ = 2 * np.pi * hinit / vinit_circ

    # we want dt to be constant for all of our orbits so we can compare
    # the velocities in the animation
    dt = tmax_circ / 180.0

    # v < v_c
    orbit_slow = earth_orbit.Trajectory(GM=G*M_E, R_crash=R_E)
    vinit_slow = 0.8*vinit_circ

    orbit_slow.integrate(vinit_slow, hinit, dt, max_rad)
    print("slow: ", orbit_slow.npts)

    # v <~ v_c
    orbit_close = earth_orbit.Trajectory(GM=G*M_E, R_crash=R_E)
    vinit_close = 0.95*vinit_circ

    orbit_close.integrate(vinit_close, hinit, dt, max_rad)
    print("close: ", orbit_close.npts)

    # v = v_c
    orbit_circ.integrate(vinit_circ, hinit, dt, max_rad)
    print("circ: ", orbit_circ.npts)

    # v_c < v < v_e
    orbit_fast= earth_orbit.Trajectory(GM=G*M_E, R_crash=R_E)
    vinit_fast = 1.2*vinit_circ

    orbit_fast.integrate(vinit_fast, hinit, dt, max_rad)
    print("fast: ", orbit_fast.npts)

    # v > v_e
    orbit_escp = earth_orbit.Trajectory(GM=G*M_E, R_crash=R_E)
    vinit_escp = 1.5*vinit_circ  # v_escp = sqrt(2) * v_circ

    orbit_escp.integrate(vinit_escp, hinit, dt, max_rad)
    print("escp: ", orbit_escp.npts)

    # plotting

    img = plt.imread("earth.png")

    # plot the orbits one by one
    iframe = 0

    trajectories = [(orbit_slow, r"$v < v_\mathrm{circular}$", "C0"),
                    (orbit_close, r"$v \lesssim v_\mathrm{circular}$", "C1"),
                    (orbit_circ, r"$v = v_\mathrm{circular}$", "C2"),
                    (orbit_fast, r"$v_\mathrm{circular} < v < v_\mathrm{escape}$", "C3"),
                    (orbit_escp, r"$v > v_\mathrm{escape}$", "C4")]

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

            max_coord = max(max(np.fabs(orbit.x[n+1]),
                                np.fabs(orbit.y[n+1])),
                            4*R_E)

            ax.axis([-max_coord, max_coord, -max_coord, 0.5 * max_coord])

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

            fig.savefig(f"escapevel_{iframe:04d}.png")

            plt.close(fig)

            iframe += 1


if __name__ == "__main__":
    orbit_plot()
