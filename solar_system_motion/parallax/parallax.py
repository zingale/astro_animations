#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt

"""
Demonstrate the principle of Parallax.  We will take Earth's orbit
to be circular.

We work in units of AU, yr, and M_sun in these units, G = 4 pi^2
"""


class ParallaxScene:
    """
    We'll treat the entire collection of the Earth/Sun, foreground
    star, and background star as an object.  The only real thing that we
    need to change from frame to frame is the location of Earth
    """

    def __init__(self):

        # start Earth on the x-axis, on the opposite side of the field of
        # stars we will reference -- we accomplish this through a phase
        self.phi = np.pi

        # number of steps per year (make this a number divisible by 4)
        self.nsteps_year = 360

        # angular velocity (radians per year)
        self.omega = 2.0*np.pi/1.0

        # semi-major axis of planet Earth
        self.a_E = 1.0

        # position of Earth over the year
        omega_t = np.arange(self.nsteps_year)*2.0*np.pi/(self.nsteps_year-1)
        self.x_orbit = self.a_E*np.cos(omega_t + self.phi)
        self.y_orbit = self.a_E*np.sin(omega_t + self.phi)

        # foreground star
        self.x_fg = 3.5
        self.y_fg = 0.0


    def draw_sun_and_orbit(self, ax=None):

        # draw the Sun
        ax.scatter([0], [0], s=1600, marker=(20, 1), color="k")
        ax.scatter([0], [0], s=1500, marker=(20, 1), color="#FFFF00")

        # plot the orbit
        ax.plot(self.x_orbit, self.y_orbit, color="C0", linestyle="--")

    def draw_earth(self, time, connect_to_fg=0, ax=None):

        x_E = self.a_E*np.cos(self.omega*time + self.phi)
        y_E = self.a_E*np.sin(self.omega*time + self.phi)

        # plot Earth
        ax.scatter([x_E], [y_E], s=100, color="C0")

        # draw the line connecting Earth and the foreground star
        if connect_to_fg == 1:
            slope = (y_E - self.y_fg)/(x_E - self.x_fg)
            xpt1 = 4.5
            ypt1 = y_E + slope*(xpt1 - x_E)
            x_E_old = x_E
            y_E_old = y_E
            ax.plot([x_E_old,xpt1], [y_E_old,ypt1], color="C2", linestyle="--")

    def draw_foreground_star(self, ax=None):
        # draw the foreground star
        ax.scatter([self.x_fg], [self.y_fg], s=200, marker=(5,1), color="C3")

        # draw the line connecting the Sun and the foreground star
        ax.plot([0,self.x_fg], [0,self.y_fg], 'k--')
        ax.text(1.5, -0.125, "d", color="k")

    def draw_background_stars(self, ax=None):
        # draw some random background stars

        pos = [(4.2, 1.6), (4.7, 1.0), (4.4, -0.4), (4.8, -0.9),
               (4.1, -1.3), (4.3, -1.5), (4.5, 0.5)]

        for x, y in pos:
            ax.scatter([x], [y], s=200, marker=(5, 1), color="C9")

    def setup_fig(self, time, name, fig=None, ax=None):

        ax.axis([-1.5, 5.0, -1.8, 1.8])

        ax.axis('off')
        ax.set_aspect("equal", "datalim")

        fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)

        fig.set_size_inches(12.8, 7.2)

        ax.set_xlabel("AU")
        ax.set_ylabel("AU")
        ax.text(-1.4, -1.8, f"time = {time:6.3f} yr")

        fig.savefig(name)


def parallax():

    t = 0.0

    nyears = 1.0

    # set the initial timestep
    nsteps_year = 512
    dt = 1.0/nsteps_year

    # compute the total number of steps needed
    nsteps = int(nyears*nsteps_year)

    p = ParallaxScene()

    iout = 0

    fig, ax = plt.subplots()

    # integrate until the Earth is at a right angle
    for n in range(nsteps//4):

        ax.clear()

        p.draw_sun_and_orbit(ax=ax)
        p.draw_earth(t, ax=ax)
        p.draw_foreground_star(ax=ax)
        p.draw_background_stars(ax=ax)

        p.setup_fig(t, f"parallax_{iout:04d}.png", fig=fig, ax=ax)

        t += dt
        iout += 1

    # show the line connecting the current position and the foreground star
    # don't advance time
    print("connecting")

    nframes = 50

    ax.clear()

    p.draw_sun_and_orbit(ax=ax)
    p.draw_earth(t, connect_to_fg=1, ax=ax)
    p.draw_foreground_star(ax=ax)
    p.draw_background_stars(ax=ax)

    plt.text(1.5, -0.8, "line of sight\nto foreground star", color="C2")

    p.setup_fig(t, f"parallax_{iout:04d}.png", fig=fig, ax=ax)
    iout += 1

    for n in range(nframes):
        fig.savefig(f"parallax_{iout:04d}.png")
        iout += 1

    # integrate for 6 months
    for n in range(nsteps//2):

        ax.clear()

        p.draw_sun_and_orbit(ax=ax)
        p.draw_earth(t, ax=ax)
        p.draw_foreground_star(ax=ax)
        p.draw_background_stars(ax=ax)

        p.setup_fig(t, f"parallax_{iout:04d}.png", fig=fig, ax=ax)

        t += dt
        iout += 1

    # show the new line connecting the current position and the foreground star
    # don't advance time
    print("connecting2")

    nframes = 50

    ax.clear()

    p.draw_sun_and_orbit(ax=ax)
    p.draw_earth(t, connect_to_fg=1, ax=ax)
    p.draw_foreground_star(ax=ax)
    p.draw_background_stars(ax=ax)

    plt.text(1.5, 1.0, "new line of sight\nto foreground star", color="C2")

    p.setup_fig(t, f"parallax_{iout:04d}.png", fig=fig, ax=ax)
    iout += 1

    for n in range(nframes):
        fig.savefig(f"parallax_{iout:04d}.png")
        iout += 1

    # integrate for the final 1/4 year
    for n in range(nsteps//4):

        ax.clear()

        p.draw_sun_and_orbit(ax=ax)
        p.draw_earth(t, ax=ax)
        p.draw_foreground_star(ax=ax)
        p.draw_background_stars(ax=ax)

        p.setup_fig(t, f"parallax_{iout:04d}.png", fig=fig, ax=ax)

        t += dt
        iout += 1

    # summarize
    nframes = nsteps_year//2
    for n in range(nframes):

        ax.clear()

        p.draw_sun_and_orbit(ax=ax)
        p.draw_earth(t, ax=ax)
        p.draw_foreground_star(ax=ax)
        p.draw_background_stars(ax=ax)


        ax.plot([0.0, 0.0], [0.0, -1.0], "C3")
        ax.plot([0.0, p.x_fg], [-1.0, p.y_fg], "C3")
        ax.text(p.x_fg-0.5, -0.1, "p", color="C2")
        ax.text(2.0, 0.5, "tan p = 1 AU / d", color="C2")
        ax.text(-0.5, -0.5, "1 AU", color="C3")

        p.setup_fig(t, f"parallax_{iout:04d}.png", fig=fig, ax=ax)

        t += dt
        iout += 1


if __name__== "__main__":
    parallax()
