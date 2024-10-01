#!/usr/bin/env python

import anim_solvers.stick_figure as sf
import anim_solvers.myarrows as arrow
import numpy as np
import matplotlib.pyplot as plt


class Earth:
    """ draw Earth and add embelishments to
        illustrate basic concepts """

    def __init__(self, latitude=40, axial_tilt=-23.5):

        self.lat = latitude

        # note that angles are measured counterclockwise from the +x
        # axis.  making it negative gives us the summer solstice
        self.tilt = axial_tilt

        # scalings for the person and earth
        self.L = 1.5
        self.R = 10

        # location of Earth
        self.x0 = 0.0
        self.y0 = 0.0

    def draw_earth(self, ax=None):

        # draw a circle
        theta = np.radians(np.arange(0,361))
        ax.plot(self.x0 + self.R*np.cos(theta),
                self.y0 + self.R*np.sin(theta), c="b")

    def draw_day_night(self, ax=None):

        theta_half = np.radians(np.arange(90,271))
        ax.fill(list(self.x0 + self.R*np.cos(theta_half)) +
                [self.x0 +self.R*np.cos(theta_half[0])],
                list(self.y0 + self.R*np.sin(theta_half)) +
                [self.y0 +self.R*np.sin(theta_half[0])], "0.75",
                zorder=-1)

    def draw_ecliptic(self, ax=None):

        L = 5.0

        ax.plot([self.x0 -1.2*self.R, self.x0 + L*self.R],
                [self.y0, self.y0], ls=":", color="k")

        ax.plot([self.x0, self.x0],
                [self.y0-self.R, self.y0+self.R], ls=":", color="k")

        ax.text(self.x0 + L*self.R, self.y0 + 0.1*self.R,
                "ecliptic", horizontalalignment="right")

    def draw_rot_axis(self, ax=None):

        overhang = 1.2

        ps = (self.x0, self.y0 - overhang*self.R)
        pe = (self.x0, self.y0 + overhang*self.R)

        ts = sf._rotate(ps, (self.x0, self.y0), np.radians(self.tilt))
        te = sf._rotate(pe, (self.x0, self.y0), np.radians(self.tilt))

        ax.plot([ts[0], te[0]], [ts[1], te[1]], color="b")

        a = arrow.ArcArrow((0, 0), 0.5*self.R,
                           theta_start=90+self.tilt, theta_end=90.0)
        a.draw(color="b", ax=ax)

        mid = 0.5*(90 + self.tilt + 90)
        ax.text(0.51*self.R*np.cos(np.radians(mid)),
                0.51*self.R*np.sin(np.radians(mid)),
                     r"$\alpha$", color="b", horizontalalignment="left")

    def draw_parallel(self, l, color="k", ls="-", label=None, ax=None):
        """ draw a line of latitude """

        # working in Earth coordinates
        # find the (x, y) of the parallel start and end
        xs = self.x0 - self.R*np.cos(np.radians(l))
        ys = self.y0 + self.R*np.sin(np.radians(l))

        xe = self.x0 + self.R*np.cos(np.radians(l))
        ye = self.y0 + self.R*np.sin(np.radians(l))

        # since we are tilted, the actual angle in the drawing
        # (from +x) is l + tilt
        angle = self.tilt

        # transform the points into the rotated frame
        ts = sf._rotate((xs, ys), (self.x0, self.y0), np.radians(angle))
        te = sf._rotate((xe, ye), (self.x0, self.y0), np.radians(angle))

        ax.plot([ts[0], te[0]], [ts[1], te[1]],
                color=color, ls=ls)

        if label is not None:
            if l > 0:
                va = "bottom"
            elif l == 0:
                va = "center"
            else:
                va = "top"

            ax.text(ts[0]-0.01*self.R, ts[1], label,
                    horizontalalignment="right",
                    verticalalignment=va, color=color)

    def draw_equator(self, ax=None):
        self.draw_parallel(0, color="b", label="equator", ax=ax)

    def draw_sun(self, ax=None):
        ax.scatter([self.x0 + 4.5*self.R], [self.y0],
                   s=2000, marker=(16, 1), zorder=100,
                   color="k")
        ax.scatter([self.x0 + 4.5*self.R], [self.y0],
                   s=1900, marker=(16, 1), zorder=100,
                   color="#FFFF00")

    def draw_tropics(self, ax=None):
        self.draw_parallel(np.abs(self.tilt), color="g", ls="--",
                           label="tropic of cancer", ax=ax)
        self.draw_parallel(-np.abs(self.tilt), color="g", ls="--",
                           label="tropic of capricorn", ax=ax)

    def draw_arctic_circles(self, ax=None):
        self.draw_parallel(90-np.abs(self.tilt), color="g",
                           ls="--", label="arctic circle",
                           ax=ax)
        self.draw_parallel(-90+np.abs(self.tilt), color="g",
                           ls="--", label="antarctic circle",
                           ax=ax)

    def draw_my_latitude(self, ax=None):

        angle = self.lat + self.tilt

        center = ((self.R + 0.5*self.L)*np.cos(np.radians(angle)),
                  (self.R + 0.5*self.L)*np.sin(np.radians(angle)))

        sf.draw_person(center, self.L, np.radians(angle - 90),
                       color="r", ax=ax)

        equator = 0 + self.tilt

        ax.plot([self.x0, self.R*np.cos(np.radians(angle))],
                [self.y0, self.R*np.sin(np.radians(angle))],
                color="r", ls="-")

        a = arrow.ArcArrow((self.x0, self.y0), 0.5*self.R,
                           theta_start=equator, theta_end=angle)
        a.draw(color="r", ax=ax)

        mid = 0.5*(equator + angle)
        ax.text(0.51*self.R*np.cos(np.radians(mid)),
                0.51*self.R*np.sin(np.radians(mid)), r"$l$", color="r",
                horizontalalignment="left")

    def draw_zenith(self, ax=None):

        angle = self.lat + self.tilt

        zenith = 3.0*self.R
        ax.plot([self.x0, zenith*np.cos(np.radians(angle))],
                [self.y0, zenith*np.sin(np.radians(angle))],
                color="r", ls=":")

        ax.text(zenith*np.cos(np.radians(angle)),
                zenith*np.sin(np.radians(angle)), "zenith", color="r",
                horizontalalignment="left")

    def draw_horizon(self, ax=None):
        angle = self.tilt + self.lat

        sc = (self.x0 + self.R*np.cos(np.radians(angle)),
              self.y0 + self.R*np.sin(np.radians(angle)))

        ps = (-0.75*self.R, 0)
        pe = (0.75*self.R, 0)

        ts = sf._rotate(ps, (0, 0), np.radians(90+angle))
        te = sf._rotate(pe, (0, 0), np.radians(90+angle))

        ax.plot([sc[0]+ts[0], sc[0]+te[0]],
                [sc[1]+ts[1], sc[1]+te[1]], color="c")
        ax.text(sc[0]+ts[0], sc[1]+ts[1], "local horizon",
                horizontalalignment="left",
                verticalalignment="top",
                rotation=270+angle, color="c")


class Scene:
    """ a container to hold the sequence of function calls we will do to
        compose the scene.  This way we can incrementally add to the
        figure """

    def __init__(self, other, xlim=None, ylim=None):
        self.other = other
        self.funcs = []
        self.xlim = xlim
        self.ylim = ylim

    def addto(self, f):
        self.funcs.append(f)

    def draw(self, description=None, ofile="test.png"):

        fig, ax = plt.subplots()
        for f in self.funcs:
            f(ax=ax)

        ax.axis("off")
        ax.set_aspect("equal", "datalim")

        fig.set_size_inches(12.8, 7.2)

        if description is not None:
            ax.text(0.025, 0.05, description,
                    transform=fig.transFigure)

        if self.xlim is not None:
            ax.set_xlim(*self.xlim)

        if self.ylim is not None:
            ax.set_ylim(*self.ylim)

        plt.subplots_adjust(left=0, bottom=0,
                            right=1, top=1, wspace=0, hspace=0)

        # dpi = 100 for 720p, 150 for 1080p
        fig.savefig(ofile, dpi=150)


def doit():

    e = Earth(latitude=42)

    sc = Scene(e, xlim=(-2*e.R, 5*e.R), ylim=(-2*e.R, 2*e.R))

    n = 0
    sc.addto(e.draw_earth)
    sc.addto(e.draw_ecliptic)
    sc.addto(e.draw_sun)
    sc.draw(ofile=f"earth_{n:02d}",
            description="Earth and the ecliptic:\n" +
            "the ecliptic is the orbital plane, connecting the Earth and the Sun")

    n += 1
    sc.addto(e.draw_day_night)
    sc.draw(ofile=f"earth_{n:02d}",
            description="the day/night line:\n" +
            "night is the hemisphere pointed away from the Sun")

    n += 1
    sc.addto(e.draw_rot_axis)
    sc.draw(ofile=f"earth_{n:02d}",
            description="Earth's axial tilt:\n" +
            r"Earth's rotation axis is tilted and angle $\alpha = 23.5^\circ$ with respect to the ecliptic")

    n += 1
    sc.addto(e.draw_equator)
    sc.draw(ofile=f"earth_{n:02d}",
            description="Earth's equator:\n" +
            "the equator is perpendicular to the rotation axis")

    n += 1
    sc.draw(ofile=f"earth_{n:02d}",
            description="the Sun:\n" +
            "the Sun is on the ecliptic (not shown to scale)\n" +
            "here Earth's North Pole is maximally pointed toward the Sun -- this is the day of the summer solstice")


    n += 1
    sc.addto(e.draw_my_latitude)
    sc.draw(ofile=f"earth_{n:02d}",
            description="latitude on Earth:\n" +
            r"latitude is just the angle above or below the equator.  Here is an observer at a latitude $l$")

    n += 1
    sc.addto(e.draw_zenith)
    sc.draw(ofile=f"earth_{n:02d}",
            description="your zenith:\n" +
            "down is the direction connecting you to the center of the Earth (the direction gravity points)\n" +
            "up is opposite down -- here the zenith is shown as the point directly above us")

    n += 1
    sc.addto(e.draw_tropics)
    sc.draw(ofile=f"earth_{n:02d}",
            description="the tropics:\n" +
            r"the tropic lines are +/- $\alpha$ in latitude -- note that the Sun is directly overhead for an observer on the Tropic of Cancer on the summer solstice")

    n += 1
    sc.addto(e.draw_arctic_circles)
    sc.draw(ofile=f"earth_{n:02d}",
            description="the arctic and antarctic circles:\n" +
            "on the summer solstice the Sun never sets between the arctic circle and North Pole -- note how everything is in daylight at these high latitudes\n" +
            "the opposite is true between the antarctic circle and the South Pole -- the Sun is never above the horizon (always night)\n" +
            r"these latitudes are just +/- $(90^\circ - \alpha)$")

    n += 1
    sc.addto(e.draw_horizon)
    sc.draw(ofile=f"earth_{n:02d}",
            description="horizon:\n" +
            "your local horizon is tangent to the surface of the Earth where you are standing")

    sc.draw()


if __name__ == "__main__":
    doit()
