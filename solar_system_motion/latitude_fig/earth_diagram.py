#!/usr/bin/env python

import anim_solvers.stick_figure as sf
import anim_solvers.myarrows as arrow
import numpy as np
import matplotlib.pylab as plt

class Earth(object):
    """ draw Earth and add embelishments to
        illustrate basic concepts """

    def __init__(self, latitude=40, axial_tilt=-23.5):

        self.lat = latitude

        # making it negative gives us the summer solstice
        self.tilt = axial_tilt

        # scalings for the person and earth
        self.L = 1.5
        self.R = 10

        # location of Earth
        self.x0 = 0.0
        self.y0 = 0.0

    def draw_earth(self, show_day_night=False):

        # draw a circle
        theta = np.radians(np.arange(0,361))
        plt.plot(self.x0 + self.R*np.cos(theta),
                 self.y0 + self.R*np.sin(theta), c="b")


        if show_day_night:
            theta_half = np.radians(np.arange(90,271))
            plt.fill(list(self.x0 + self.R*np.cos(theta_half)) +
                     [self.x0 +self.R*np.cos(theta_half[0])],
                     list(self.y0 + self.R*np.sin(theta_half)) +
                     [self.y0 +self.R*np.sin(theta_half[0])], "0.75",
                     zorder=-1)


    def draw_ecliptic(self):

        L = 5.0

        plt.plot([self.x0 -1.2*self.R, self.x0 + L*self.R],
                 [self.y0, self.y0], ls=":", color="k")

        plt.plot([self.x0, self.x0],
                 [self.y0-self.R, self.y0+self.R], ls=":", color="k")

        plt.text(self.x0 + L*self.R, self.y0 + 0.1*self.R,
                 "ecliptic", horizontalalignment="right")


    def draw_rot_axis(self):

        overhang = 1.1

        ps = (self.x0, self.y0 - overhang*self.R)
        pe = (self.x0, self.y0 + overhang*self.R)

        ts = sf._rotate(ps, (self.x0, self.y0), np.radians(self.tilt))
        te = sf._rotate(pe, (self.x0, self.y0), np.radians(self.tilt))

        plt.plot([ts[0], te[0]], [ts[1], te[1]], color="b")


        a = arrow.ArcArrow((0, 0), 0.5*self.R, theta_start=90+self.tilt, theta_end=90.0)
        a.draw(color="b")

        mid = 0.5*(90 + self.tilt + 90)
        plt.text(0.51*self.R*np.cos(np.radians(mid)),
                 0.51*self.R*np.sin(np.radians(mid)), r"$\alpha$", color="b", horizontalalignment="left")


    def draw_equator(self):

        ps = (self.x0, self.y0 - self.R)
        pe = (self.x0, self.y0 + self.R)

        ts = sf._rotate(ps, (self.x0, self.y0),
                        np.radians(90+self.tilt))
        te = sf._rotate(pe, (self.x0, self.y0),
                        np.radians(90+self.tilt))

        plt.plot([ts[0], te[0]], [ts[1], te[1]], color="b")


    def draw_latitude(self):

        angle = self.lat + self.tilt

        center = ( (self.R + 0.5*self.L)*np.cos(np.radians(angle)),
                   (self.R + 0.5*self.L)*np.sin(np.radians(angle)) )
        sf.draw_person(center, self.L, np.radians(angle - 90), color="r")


        zenith = 3.0*self.R
        plt.plot([0.0, zenith*np.cos(np.radians(angle))],
                 [0.0, zenith*np.sin(np.radians(angle))], color="r", ls=":")
        plt.text(zenith*np.cos(np.radians(angle)),
                 zenith*np.sin(np.radians(angle)), "zenith", color="r",
                 horizontalalignment="left")

        equator = 0 + self.tilt

        a = arrow.ArcArrow((0, 0), 0.5*self.R, theta_start=equator, theta_end=angle)
        a.draw(color="r")

        mid = 0.5*(equator + angle)
        plt.text(0.51*self.R*np.cos(np.radians(mid)),
                 0.51*self.R*np.sin(np.radians(mid)), r"$l$", color="r", horizontalalignment="left")

    def draw_zenith(self):
        pass



def doit():

    e = Earth(latitude=42)

    e.draw_earth(show_day_night=True)
    e.draw_ecliptic()
    e.draw_rot_axis()
    e.draw_equator()
    e.draw_latitude()

    plt.axis("off")

    ax = plt.gca()
    ax.set_aspect("equal", "datalim")

    f = plt.gcf()
    f.set_size_inches(12.8, 7.2)

    plt.savefig("test.png")



if __name__ == "__main__":
    doit()
