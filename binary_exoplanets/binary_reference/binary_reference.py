#!/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

import anim_solvers.binary_integrator as bi

# compute the orbits of stars in a binary system.
#
# Here we put the star 1 at the origin -- this shows what the orbit of
# the second star looks like from the point of view of the first.
#
# This version allows for elliptical orbits with some arbitrary
# orientation wrt to the observer (although, still face-on)
#
# The plotting assumes that M_2 < M_1
#
# M. Zingale (2009-02-12)

# we work in CGS units

M_sun = bi.M_sun
AU = bi.AU
year = bi.year
G = bi.G


# a simple class to serve as a container for the orbital information
# for the two stars

def find_scinotat(number):

    b = int(np.log10(np.fabs(number)))
    a = number/10**b

    return a, b


def radial_velocity():

    # set the masses
    M_star1 = 4.0*M_sun      # star 1's mass
    M_star2 = M_sun      # star 2's mass

    # set the semi-major axis of the star 2 (and derive that of star 1)
    # M_star2 a_star2 = -M_star1 a_star1 (center of mass)
    a_star2 = 10.0*AU
    a_star1 = (M_star2/M_star1)*a_star2

    # set the eccentricity
    ecc = 0.4

    # set the angle to rotate the semi-major axis wrt the observer
    theta = np.pi/6.0

    # display additional information
    annotate = True

    # create the binary object
    b = bi.Binary(M_star1, M_star2, a_star1 + a_star2, ecc, theta,
                  annotate=annotate)

    # set the timestep in terms of the orbital period
    dt = b.P/360.0
    tmax = 2.0*b.P  # maximum integration time

    b.integrate(dt, tmax)
    s1 = b.orbit1
    s2 = b.orbit2

    # ================================================================
    # plotting
    # ================================================================

    iframe = 0

    xmin = 1.05*(s2.x - s1.x).min()
    xmax = 1.05*(s2.x - s1.x).max()

    ymin = 1.05*(s2.y - s1.y).min()
    ymax = 1.05*(s2.y - s1.y).max()

    for n in range(len(s1.t)):

        fig = plt.figure(1)
        fig.clear()

        plt.subplots_adjust(left=0.025, right=0.975, bottom=0.025, top=0.95)

        axl = fig.add_subplot(121)
        axr = fig.add_subplot(122)

        for i in range(2):

            if i == 0:
                ax = axl

                # offsets -- we want star 1 at the origin
                xoffset = s1.x[n]
                yoffset = s1.y[n]

            else:
                ax = axr

                # no offsets
                xoffset = 0.0
                yoffset = 0.0

            ax.set_aspect("equal", "datalim")
            ax.set_axis_off()

            # center of mass
            ax.scatter([0-xoffset], [0-yoffset], s=150, marker="x", color="k")

            if i == 0:
                # plot star 1's position
                symsize = 200
                ax.scatter([s1.x[n]-xoffset], [s1.y[n]-yoffset], s=symsize, color="C0")

                # plot star 2's orbit and position
                symsize = 200*(b.M2/b.M1)
                ax.plot(s2.x-s1.x, s2.y-s1.y, color="C1")

                ax.scatter([s2.x[n]-xoffset], [s2.y[n]-yoffset], s=symsize, color="C1", zorder=100)

            else:
                # plot star 1's orbit position
                symsize = 200
                ax.scatter([s1.x[n]], [s1.y[n]], s=symsize, color="C0", zorder=100)
                ax.plot(s1.x, s1.y, color="C0")

                # plot star 2's orbit and position
                symsize = 200*(b.M2/b.M1)
                ax.scatter([s2.x[n]], [s2.y[n]], s=symsize, color="C1", zorder=100)
                ax.plot(s2.x, s2.y, color="C1")

            if i == 0 and annotate:
                # display time
                ax.text(0.05, 0.05, "time = {:6.3f} yr".format(s1.t[n]/year),
                        transform=ax.transAxes)

                # display information about stars
                ax.text(0.025, 0.9, r"mass ratio: {:3.2f}".format(b.M1/b.M2),
                        transform=ax.transAxes, color="k", fontsize="medium")
                ax.text(0.025, 0.86, r"eccentricity: {:3.2f}".format(b.e),
                        transform=ax.transAxes, color="k", fontsize="medium")



            if annotate:
                # plot a reference line
                ax.plot([0, 1*AU], [0.93*ymin, 0.93*ymin], color="k")
                ax.text(0.5*AU, 0.975*ymin, "1 AU",
                        horizontalalignment="center", verticalalignment="top")

            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)

            if i == 0:
                ax.set_title("massive star frame of reference")
            else:
                ax.set_title("center of mass frame of reference")

        fig.set_size_inches(12.8, 7.2)

        plt.savefig("binary_star_{:04d}.png".format(iframe))

        iframe += 1


if __name__== "__main__":
    radial_velocity()
