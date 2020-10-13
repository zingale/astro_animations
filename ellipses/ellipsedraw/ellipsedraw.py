#!/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axisartist.axislines import SubplotZero

# show how ellipses are drawn

# M. Zingale

def ellipse():

    # theta ranges from 0 to 2pi
    npts = 360
    theta = np.linspace(0, 2*np.pi, npts)

    e = 0.5
    a = 1.0

    r = a*(1.0 - e**3)/(1.0 + e*np.cos(theta))

    x = r*np.cos(theta)
    y = r*np.sin(theta)


    # plotting
    for n in range(npts):

        fig = plt.figure(1)
        fig.clear()

        ax = SubplotZero(fig, 111)
        fig.add_subplot(ax)

        ax.set_aspect("equal", "datalim")

        ax = plt.gca()
        ax.set_aspect("equal", "datalim")

        # draw the ellipse
        ax.plot(x, y, color="k", linewidth=2)

        # draw our current point
        ax.scatter([x[n]], [y[n]], color="k", s=75)

        # second foci
        ax.scatter([-2.0*a*e], [0], color="C0", marker="x", s=200)

        # primary foci
        ax.scatter([0], [0], color="C1", marker="x", s=200)

        # draw lines connecting the foci to the current point
        ax.plot([0, x[n]], [0, y[n]], color="C1", zorder=100, linewidth=2)
        ax.plot([-2.0*a*e, x[n]], [0, y[n]], color="C0", zorder=100, linewidth=2)

        ax.set_xlim(-2.5, 1.5)
        ax.set_ylim(-2., 2.)

        len1 = np.sqrt((x[n] - 0)**2 + (y[n] - 0)**2)
        len2 = np.sqrt((x[n] - (-2.0*a*e))**2 + (y[n] - 0)**2)

        ax.set_title("Ellipse, eccenticity = {:5.3f}".format(e))

        ax.text(-1.5, -1.25, "r length: {:5.3f}".format(len1), color="C1")
        ax.text(-1.5, -1.5, "r' length: {:5.3f}".format(len2), color="C0")
        ax.text(-1.5, -1.75, "r + r' = {:5.3f}".format(len1 + len2))

        for direction in ["xzero", "yzero"]:
            # adds arrows at the ends of each axis
            ax.axis[direction].set_axisline_style("-|>")

            # adds X and Y-axis from the origin
            ax.axis[direction].set_visible(True)

        for direction in ["left", "right", "bottom", "top"]:
            # hides borders
            ax.axis[direction].set_visible(False)

        fig.set_size_inches(7.2, 7.2)

        fig.set_size_inches(7.2,7.2)

        plt.savefig("ellipsedraw_{:03d}.png".format(n))

if __name__== "__main__":
    ellipse()
