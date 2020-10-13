#!/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

import anim_solvers.binary_integrator as bi

# compute the orbits of stars in a binary system.
#
# Here we put the center of mass at the origin.
#
# This version allows for elliptical orbits with some arbitrary
# orientation wrt to the observer (although, still face-on)
#
# The plotting assumes that M_2 < M_1
#
# M. Zingale (2009-02-12)

# we work in CGS units
G = bi.G
M_sun = bi.M_sun
AU = bi.AU
year = bi.year


def astrometric_binary():

    # set the masses
    M_star1 = M_sun           # star 1's mass
    M_star2 = 0.5*M_sun      # star 2's mass

    # set the semi-major axis of the star 2 (and derive that of star 1)
    # M_star2 a_star2 = -M_star1 a_star1 (center of mass)
    a_star2 = 1.0*AU
    a_star1 = (M_star2/M_star1)*a_star2

    # set the eccentricity
    ecc = 0.0

    # set the angle to rotate the semi-major axis wrt the observer
    theta = np.pi/2.0

    b = bi.Binary(M_star1, M_star2, a_star1+a_star2, ecc, theta)

    # velocity of the center of mass
    v = 10*a_star2/b.P


    # set the timestep in terms of the orbital period
    dt = b.P/360.0
    tmax = 2.0*b.P  # maximum integration time

    b.integrate(dt, tmax)
    s1 = b.orbit1
    s2 = b.orbit2

    # now move the center of mass according to the velocity defined
    # above, and update the positions
    x_cm = np.zeros(len(s1.t), np.float64)
    y_cm = np.zeros(len(s2.t), np.float64)

    for n in range(len(s1.t)):

        x_cm[n] = v*s1.t[n]
        y_cm[n] = 0.0

        s1.x[n] += v*s1.t[n]
        s2.x[n] += v*s2.t[n]


    # ================================================================
    # plotting
    # ================================================================

    fig = plt.figure(1)
    fig.clear()

    plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)

    ax = fig.add_subplot(111)
    ax.set_aspect("equal", "datalim")
    ax.set_axis_off()

    # plot the center of mass motion
    ax.scatter([0], [0], s=50, marker="x", color="k")
    ax.plot(x_cm, y_cm, color="0.5", linestyle="--")

    # plot star 1's orbit and position
    ax.plot(s1.x, s1.y, color="C1")
    ax.scatter([s1.x[0]],[s1.y[0]], s=200, color="C1")

    # plot star 2's orbit and position
    ax.plot(s2.x, s2.y, color="C0")
    ax.scatter([s2.x[0]], [s2.y[0]], s=100, color="C0")

    ax.text(10*a_star2, 3*a_star2, "mass ratio: %3.2f" % (b.M1/b.M2),
            color="k", verticalalignment="center")
    ax.text(10*a_star2, 2.5*a_star2, "eccentricity: %3.2f" % (ecc),
            color="k", verticalalignment="center")

    # labels
    ax.plot([0.1*a_star2,1.1*a_star2], [2.9*a_star2,2.9*a_star2],
            color="0.5", linestyle="--")
    ax.text(1.3*a_star2, 2.9*a_star2, "path of center of mass",
            verticalalignment="center", color="0.5")

    ax.plot([0.1*a_star2,1.1*a_star2], [2.4*a_star2,2.4*a_star2], color="C1")
    ax.text(1.3*a_star2, 2.4*a_star2, "path of primary star",
            verticalalignment="center", color="C1")

    ax.plot([0.1*a_star2,1.1*a_star2],[1.9*a_star2,1.9*a_star2], color="C0")
    ax.text(1.3*a_star2, 1.9*a_star2, "path of unseen companion",
            verticalalignment="center", color="C0")

    ax.set_xlim(-0.3*a_star2, 19.7*a_star2)
    ax.set_ylim(-3*a_star2, 3*a_star2)

    fig.set_size_inches(10.0,3.0)

    plt.savefig("astrometric_binary.png")


if __name__== "__main__":
    astrometric_binary()
