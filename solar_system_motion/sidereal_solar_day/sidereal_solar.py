import numpy as np
import matplotlib.pyplot as plt
import anim_solvers.solar_system_integrator as ssi
from anim_solvers.stick_figure import draw_person

"""
Illustrate the difference between the sidereal and solar day for earth
we exaggerate the effect to make it clearer
"""


def sidereal():

    # set the semi-major axis and eccentricity
    a = 1.0
    e = 0.0

    ss = ssi.SolarSystem()

    x = 0.0
    y = -a*(1.0 + e)
    vx = np.sqrt((ss.GM/a)*(1.0+e)/(1.0-e))
    vy = 0.0

    pos = ssi.PlanetPosVel(x, y, vx, vy)

    ss.add_planet(a, e, loc="specify", pos_vel=pos)

    # compute the period of the orbit from Kepler's law
    P_orbital = ss.period(0)

    # set the rotation period -- this is the sidereal day
    #P_rotation = 0.99726968*day    # Earth
    P_rotation = 0.05*P_orbital

    omega = 2*np.pi/P_rotation

    # compute the length of the solar day
    P_solar = P_rotation/(1.0 - P_rotation/P_orbital)

    print("period = ", P_orbital)

    print("sidereal day = ", P_rotation)
    print("solar day    = ", P_solar)

    # We will maintain two trajectories: full_orbit is the full orbit
    # of the planet.  orbit is just a small segment carrying 2
    # sidereal days

    # evolve for 2 sidereral days with lots of frames to get a smooth
    # animation
    tmax = 2.0*P_rotation
    nsteps = 1000

    # fraction of a year we are doing
    f = tmax/P_orbital
    orbit = ss.integrate(nsteps/f, f)
    dt = orbit[0].t[1] - orbit[0].t[0]

    # evolve the full orbit for plotting purposes
    nsteps_per_year = 720.0
    full_orbit = ss.integrate(nsteps_per_year, 1.0)

    # plotting

    # plot the orbit
    iframe = 0
    for n in range(len(orbit[0].t)):

        fig, ax = plt.subplots()

        # plot the foci
        ax.scatter([0], [0], s=1600, marker=(20, 1), color="k")
        ax.scatter([0], [0], s=1500, marker=(20, 1), color="#FFFF00")

        # plot the orbit
        ax.plot(full_orbit[0].x, full_orbit[0].y, color="0.5")

        # plot planet
        theta = np.radians(np.arange(360))
        r = 0.075  # exaggerate the planet's size
        x_surface = orbit[0].x[n] + r*np.cos(theta)
        y_surface = orbit[0].y[n] + r*np.sin(theta)
        ax.fill(x_surface, y_surface, "C0", edgecolor="C0", zorder=1000)

        # plot a point on the planet's surface
        xpt = orbit[0].x[n] + 1.45*r*np.cos(omega*orbit[0].t[n]+np.pi/2.0)
        ypt = orbit[0].y[n] + 1.45*r*np.sin(omega*orbit[0].t[n]+np.pi/2.0)
        #ax.scatter([xpt], [ypt], s=15, color="k")
        draw_person([xpt, ypt], 0.075, omega*orbit[0].t[n], ax=ax)

        ax.axis([-1.15, 1.15, -1.15, 1.15])
        ax.axis("off")
        ax.set_aspect("equal", "datalim")

        fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)
        fig.set_size_inches(7.2, 7.2)

        ax.set_title("Sidereal vs. Solar Day")

        ax.text(-1.0, -1.2, "Note: length of day exaggerated",
                fontsize=9, color="0.5")

        if (orbit[0].t[n] >= P_rotation - 0.5*dt and
            orbit[0].t[n] <= P_rotation + 0.5*dt):
            n_sidereal = n

        # special lines
        if (orbit[0].t[n] >= P_rotation - 0.5*dt and
            orbit[0].t[n] <= P_solar + 0.5*dt):
            ax.plot([orbit[0].x[n_sidereal], orbit[0].x[n_sidereal]],
                    [orbit[0].y[n_sidereal], 0], linestyle=":", color="0.5")

        # special notes
        if (orbit[0].t[n] < 1.5*dt):
            ax.text(-1.0, -1.15, "Noon (Sun is on the meridian)")

            for d in range(150):
                fig.savefig(f"earth_{iframe:04d}.png")
                iframe += 1

        elif (orbit[0].t[n] >= P_rotation - 0.5*dt and
              orbit[0].t[n] <= P_rotation + 0.5*dt):
            ax.text(-1.0, -1.15, "1 Sidereal period (Earth rotated 360 degrees)")

            for d in range(150):
                fig.savefig(f"earth_{iframe:04d}.png")
                iframe += 1

        elif (orbit[0].t[n] >= P_solar - 0.5*dt and
              orbit[0].t[n] <= P_solar + 0.5*dt):
            ax.text(-1.0, -1.15, "1 Solar period (Sun is back on the meridian)")

            for d in range(150):
                fig.savefig(f"earth_{iframe:04d}.png")
                iframe += 1

        fig.savefig(f"earth_{iframe:04d}.png")
        plt.close(fig)
        iframe += 1


if __name__ == "__main__":
    sidereal()
