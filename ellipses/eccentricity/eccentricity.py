#!/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axisartist.axislines import SubplotZero

# draw ellipses of varying eccentricity

# M. Zingale

def ellipse():

    # theta ranges from 0 to 2pi
    npts = 360
    theta = np.linspace(0, 2*np.pi, npts)

    # well go through a range of eccentricities (forwards and backwards)
    n_ecc = 200
    e_tmp = np.linspace(0, 0.95, n_ecc)

    ecc = np.zeros(2*n_ecc)
    ecc[0:n_ecc] = e_tmp[:]
    ecc[n_ecc:] = e_tmp[::-1]

    a = 1.0

    for n, e in enumerate(ecc):

        r = a*(1.0 - e**2)/(1.0 + e*np.cos(theta))

        x = r*np.cos(theta)
        y = r*np.sin(theta)

        # plotting
        fig = plt.figure(1)
        fig.clear()

        ax = SubplotZero(fig, 111)
        fig.add_subplot(ax)

        ax.set_aspect("equal", "datalim")

        ax.plot(x, y, color="k", linewidth=2)

        # second foci
        ax.scatter([-2.0*a*e], [0], color="C0", marker="x", s=100)

        # primary foci
        ax.scatter([0], [0],color="C1", marker="x", s=100)

        ax.set_xlim(-2.5, 1.5)
        ax.set_ylim(-2., 2.)

        ax.text(-2.0, -1.5, f"a = {a:5.3f}, e = {e:6.4f}")

        for direction in ["xzero", "yzero"]:
            # adds arrows at the ends of each axis
            ax.axis[direction].set_axisline_style("-|>")

            # adds X and Y-axis from the origin
            ax.axis[direction].set_visible(True)

        for direction in ["left", "right", "bottom", "top"]:
            # hides borders
            ax.axis[direction].set_visible(False)

        fig.set_size_inches(7.2, 7.2)

        plt.savefig(f"ellipse_{n:03d}.png")

if __name__== "__main__":
    ellipse()
