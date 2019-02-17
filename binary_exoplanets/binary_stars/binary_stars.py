#!/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt

import anim_solvers.binary_integrator as bi


M_sun = bi.M_sun
AU = bi.AU
year = bi.year
G = bi.G

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

def find_scinotat(number):

    b = int(np.log10(np.fabs(number)))
    a = number/10**b

    return a, b


def radial_velocity(M1=1, M2=1, a2=10, e=0.0, annotate=False):

    # set the masses
    M_star1 = M1*M_sun      # star 1's mass
    M_star2 = M2*M_sun      # star 2's mass

    # set the semi-major axis of the star 2 (and derive that of star 1)
    # M_star2 a_star2 = -M_star1 a_star1 (center of mass)
    a_star2 = a2*AU
    a_star1 = (M_star2/M_star1)*a_star2

    # set the eccentricity
    ecc = e

    # set the angle to rotate the semi-major axis wrt the observer
    theta = np.pi/6.0

    # create the binary object
    b = bi.Binary(M_star1, M_star2, a_star1 + a_star2, ecc, theta, annotate=annotate)


    # set the timestep in terms of the orbital period
    dt = b.P/360.0
    tmax = 2.0*b.P  # maximum integration time

    s1, s2 = b.integrate(dt, tmax)


    # ================================================================
    # plotting
    # ================================================================

    iframe = 0

    for n in range(len(s1.t)):

        fig = plt.figure(1)
        fig.clear()

        ax = fig.add_subplot(111)

        plt.subplots_adjust(left=0.025, right=0.975, bottom=0.025, top=0.975)

        ax.set_aspect("equal", "datalim")
        ax.set_axis_off()

        ax.scatter([0], [0], s=150, marker="x", color="k")

        # if e = 0 and M_star1 = M_star2, then the orbits lie on top of one
        # another, so plot only a single orbital line.

        # plot star 1's orbit and position
        symsize = 200
        if not (b.M1 == b.M2 and b.e == 0.0):
            ax.plot(s1.x, s1.y, color="C0")
        else:
            ax.plot(s1.x, s1.y, color="k")

        ax.scatter([s1.x[n]], [s1.y[n]], s=symsize, color="C0", zorder=100)

        # plot star 2's orbit and position
        symsize = 200*(b.M2/b.M1)
        if not (b.M1 == b.M2 and b.e == 0.0):
            ax.plot(s2.x, s2.y, color="C1")

        ax.scatter([s2.x[n]], [s2.y[n]], s=symsize, color="C1", zorder=100)

        if annotate:
            # plot a reference line
            ax.plot([0, 1*AU], [-1.2*b.a2, -1.2*b.a2], color="k")
            ax.text(0.5*AU, -1.4*b.a2, "1 AU",
                    horizontalalignment='center')

            # display time
            ax.text(0.05, 0.05, "time = {:6.3f} yr".format(s1.t[n]/year),
                    transform=ax.transAxes)

        # display information about stars
        ax.text(0.05, 0.95, r"mass ratio: {:3.2f}".format(b.M1/b.M2),
                transform=ax.transAxes, color="k", fontsize="large")
        ax.text(0.05, 0.9, r"eccentricity: {:3.2f}".format(b.e),
                transform=ax.transAxes, color="k", fontsize="large")

        # energies
        if annotate:
            KE1, KE2 = b.kinetic_energies()
            PE = b.potential_energy()

            print(KE1, KE2, PE, KE1 + KE2 + PE)

            # KE 1
            sig, ex = find_scinotat(KE1)
            ax.text(0.05, 0.4, r"$K_1 = {:+4.2f} \times 10^{{{:2d}}}$ erg".format(sig, ex),
                    transform=ax.transAxes, color="C0")

            sig, ex = find_scinotat(KE2)
            ax.text(0.05, 0.35, r"$K_2 = {:+4.2f} \times 10^{{{:2d}}}$ erg".format(sig, ex),
                    transform=ax.transAxes, color="C1")

            sig, ex = find_scinotat(PE)
            ax.text(0.05, 0.3, r"$U = {:+4.2f} \times 10^{{{:2d}}}$ erg".format(sig, ex),
                    transform=ax.transAxes)

            sig, ex = find_scinotat(KE1 + KE2 + PE)
            ax.text(0.05, 0.25, r"$E = {:+4.2f} \times 10^{{{:2d}}}$ erg".format(sig, ex),
                    transform=ax.transAxes)

        #ax.set_xlim(-2.0*b.a2, 2.0*b.a2)
        #ax.set_ylim(-1.75*b.a2, 1.75*b.a2)

        fig.set_size_inches(12.8, 7.2)

        plt.tight_layout()
        plt.savefig("binary_star_mratio={:3.2f}_e={:3.2f}_{:04d}.png".format(M_star1/M_star2, e, iframe))

        iframe += 1


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-a", help="the semi-major axis of star 2", type=float, default=10.0)
    parser.add_argument("-e", help="the eccentricity", type=float, default=0.0)
    parser.add_argument("--mass1", help="mass of star 1", type=float, default=1.0)
    parser.add_argument("--mass2", help="mass of star 2", type=float, default=1.0)
    parser.add_argument("--annotate", help="--show energy details", action="store_true")

    args = parser.parse_args()

    radial_velocity(M1=args.mass1, M2=args.mass2, e=args.e, a2=args.a, annotate=args.annotate)
