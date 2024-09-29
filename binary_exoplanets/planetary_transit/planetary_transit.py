"""
compute the orbits of a transiting planet, seen nearly edge-on.

Here we put the center of mass at the origin.

This version allows for elliptical orbits with some arbitrary
orientation wrt to the observer (although, still face-on)

The plotting assumes that M_p (2) << M_star (1)
"""


import numpy as np
import matplotlib.pyplot as plt

from anim_solvers.binary_integrator import Binary

# we work in CGS units
G = 6.67428e-8        # cm^3 g^{-1} s^{-2}
M_sun = 1.98892e33    # g
AU = 1.49598e13       # cm
year = 3.1556926e7


def radial_velocity():

    # set the masses
    M_star1 = M_sun         # star 1's mass
    M_star2 = 0.2*M_sun     # planet's mass (exaggerated)

    # set the semi-major axis of the star 2 (and derive that of star 1)
    # M_star2 a_star2 = -M_star1 a_star1 (center of mass)
    a_star2 = 0.5*AU
    a_star1 = (M_star2/M_star1)*a_star2
    a = a_star1 + a_star2

    # set the eccentricity
    ecc = 0.0

    # set the angle to rotate the semi-major axis wrt the observer
    theta = 0.0

    # set the incliination wrt the observer.
    inc = 88.75    # degrees

    # create the solar system container
    ss = Binary(M_star1, M_star2, a, ecc, theta)

    # set the radii and temperatures -- dimensionless -- we are going
    # to exaggerate their scale
    R1 = 20
    R2 = 2

    rad_scal = 0.01*AU

    # set the temperatures
    T1 = 5000 # K
    T2 = 1000  # K

    # set the timestep in terms of the orbital period
    dt = ss.P / 360.0
    tmax = ss.P  # maximum integration time

    ss.integrate(dt, tmax)

    # apply the projection to account for the inclination wrt the
    # observer
    ss.orbit1.y[0:ss.npts] *= np.cos(np.radians(inc))
    ss.orbit2.y[0:ss.npts] *= np.cos(np.radians(inc))

    # we will keep one star fixed and draw the relative `orbit' of the
    # other
    x_rel = ss.orbit1.x - ss.orbit2.x
    y_rel = ss.orbit1.y - ss.orbit2.y

    # cut angle -- we don't want to plot the orbit where the stationary
    # star is (origin), so specify the angle wrt the +y axis where we
    # don't plot
    cut_angle = np.radians(20)

    frac = cut_angle/(2.0*np.pi)

    # star 2 starts on the -x axis.  3/4 through the orbit, it will
    # be behind the star 1.  Compute the range of steps to skip plotting
    cut_index1 = int(0.75 * ss.npts - frac * ss.npts)
    cut_index2 = int(0.75 * ss.npts + frac * ss.npts)

    # light curve

    # first compute the fluxes -- since we are going to normalize,
    # don't worry about the pi and sigma (Stefan-Boltzmann constant)

    # here we are assuming that R2 < R1
    f1 = R1**2 * T1**4
    f2 = R2**2 * T2**4

    f_normal = f1 + f2
    f_star2_transit = (R1**2 - R2**2) * T1**4 + R2**2 * T2**4
    f_star2_blocked = f1

    # relative fluxes -- normalized to normal
    f_star2_transit = f_star2_transit / f_normal
    f_star2_blocked = f_star2_blocked / f_normal

    print("f_star2_transit = ", f_star2_transit)
    print("f_star2_blocked = ", f_star2_blocked)

    # determine the times of the eclipses / transits
    # t_a = star 2 begins to pass in front of star 1
    # t_b = star 2 fully in front of star 1
    # t_c = star 2 begins to finish its transit of star 1
    # t_d = star 2 fully finished with transit
    # t_e = star 2 begins to go behind star 1
    # t_f = star 2 fully behind star 1
    # t_g = star 2 begins to emerge from behind star 1
    # t_h = star 2 fully emerged from behind star 1
    t_a = t_b = t_c = t_d = t_e = t_f = t_g = t_h = -1.0

    for n in range(ss.npts):

        if y_rel[n] <= 0:

            # star 2 in front of star 1
            if x_rel[n] + R2*rad_scal > -R1*rad_scal and t_a == -1:
                t_a = ss.orbit1.t[n]

            if x_rel[n] - R2*rad_scal >= -R1*rad_scal and t_b == -1:
                t_b = ss.orbit1.t[n]

            if x_rel[n] + R2*rad_scal > R1*rad_scal and t_c == -1:
                t_c = ss.orbit1.t[n]

            if x_rel[n] - R2*rad_scal >= R1*rad_scal and t_d == -1:
                t_d = ss.orbit1.t[n]

        else:

            # star 2 behind star 1
            if x_rel[n] - R2*rad_scal < R1*rad_scal and t_e == -1:
                t_e = ss.orbit1.t[n]

            if x_rel[n] + R2*rad_scal <= R1*rad_scal and t_f == -1:
                t_f = ss.orbit1.t[n]

            if x_rel[n] - R2*rad_scal < -R1*rad_scal and t_g == -1:
                t_g = ss.orbit1.t[n]

            if x_rel[n] + R2*rad_scal <= -R1*rad_scal and t_h == -1:
                t_h = ss.orbit1.t[n]

    # make an array of the flux vs. time -- this is the light curve
    f_system = np.zeros(ss.npts, np.float64)

    for n in range(ss.npts):

        f_system[n] = 1.0

        # star 2 passing in front of star 1

        if ss.orbit1.t[n] >= t_a and ss.orbit1.t[n] < t_b:
            # linearly interpolate between f = 1 at t = t_a and
            # f = f_star2_transit at t = t_b
            slope = (f_star2_transit - 1.0)/(t_b - t_a)
            f_system[n] = slope*(ss.orbit1.t[n] - t_b) + f_star2_transit

        elif ss.orbit1.t[n] >= t_b and ss.orbit1.t[n] < t_c:
            f_system[n] = f_star2_transit

        elif ss.orbit1.t[n] >= t_c and ss.orbit1.t[n] < t_d:
            # linearly interpolate between f = f_star2_transit at
            # t = t_c and f = 1 at t = t_d
            slope = (1.0 - f_star2_transit)/(t_d - t_c)
            f_system[n] = slope*(ss.orbit1.t[n] - t_d) + 1.0

        # star 2 passing behind star 1

        elif ss.orbit1.t[n] >= t_e and ss.orbit1.t[n] < t_f:
            # linearly interpolate between f = 1 at t = t_e and
            # f = f_star2_blocked at t = t_f
            slope = (f_star2_blocked - 1.0)/(t_f - t_e)
            f_system[n] = slope*(ss.orbit1.t[n] - t_f) + f_star2_blocked

        elif ss.orbit1.t[n] >= t_f and ss.orbit1.t[n] < t_g:
            f_system[n] = f_star2_blocked

        elif ss.orbit1.t[n] >= t_g and ss.orbit1.t[n] < t_h:
            # linearly interpolate between f = f_star2_blocked
            # at t = t_g and f = 1 at t = t_h
            slope = (1.0 - f_star2_blocked)/(t_h - t_g)
            f_system[n] = slope*(ss.orbit1.t[n] - t_h) + 1.0

    # plotting

    iframe = 0

    for n in range(ss.npts):

        fig = plt.figure()
        plt.subplots_adjust(left=0.15, right=0.9, bottom=0.1, top=0.9)

        ax = fig.add_subplot(211)
        ax.set_aspect("equal", "datalim")
        ax.axis("off")

        if y_rel[n] < 0.0:
            # star 2 is in front of star 1

            # plot star 1's orbit position -- set to be the origin
            xc1, yc1 = circle(0, 0, R1*rad_scal)
            ax.fill(xc1, yc1, 'C1', ec="none")

            # plot star 2's relative orbit and position
            xc2, yc2 = circle(x_rel[n], y_rel[n], R2*rad_scal)
            ax.fill(xc2, yc2, 'C0', ec="none")

        else:
            # star 1 is in front of star 2

            # plot star 2's relative orbit and position
            xc2, yc2 = circle(x_rel[n], y_rel[n], R2*rad_scal)
            ax.fill(xc2, yc2, 'C0', ec="none")

            # plot star 1's orbit position -- set to be the origin
            xc1, yc1 = circle(0, 0, R1*rad_scal)
            ax.fill(xc1, yc1, 'C1', ec="none")

        # plot the orbit -- in two segment
        ax.plot(x_rel[0:cut_index1], y_rel[0:cut_index1],
                color="0.5", linestyle="--", alpha=0.5)
        ax.plot(x_rel[cut_index2:ss.npts], y_rel[cut_index2:ss.npts],
                color="0.5", linestyle="--", alpha=0.5)

        ax.axis([-1.5*a_star2, 1.5*a_star2, -1.25*a_star2, 1.25*a_star2])

        ax = fig.add_subplot(212)

        ax.plot(ss.orbit1.t[0:ss.npts]/tmax, f_system, "k")
        ax.scatter([ss.orbit1.t[n]/tmax], [f_system[n]], color="k")

        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.9, 1.1)

        ax.set_xlabel("t/P")
        ax.set_ylabel("relative flux")

        fig.set_size_inches(6.0, 7.2)

        outfile = f"planet_transit_{iframe:04d}.png"
        fig.savefig(outfile)
        plt.close(fig)

        iframe += 1


def circle(xc, yc, r):

    theta = np.arange(361)*2.0*np.pi/360

    x = r*np.cos(theta) + xc
    y = r*np.sin(theta) + yc

    return x, y


if __name__ == "__main__":
    radial_velocity()
