#!/bin/env python

import math
import numpy as np
import matplotlib.pyplot as plt

# Show how the time between full moons is longer than the sidereal orbital
# period of the moon because of the motion of the Earth around the Sun.

# M. Zingale (2010-12-16)

# periods of the bodies
s_per_day = 24*3600.         # seconds in a day
P_E = 365.25*s_per_day       # orbital period of Earth
P_M_sid = 27.32*s_per_day    # sidereal period of Moon
P_M_syn = 29.53*s_per_day    # synodic period of Moon

# compute angular velocities
omega_E = 2.0*math.pi/P_E
omega_M = 2.0*math.pi/P_M_sid

# orbital radii (dimensionless and scale exaggerated)
r_E = 1.0
r_M = 0.25

# center of Earth
x_E_c = 0.0
y_E_c = 0.0

# starting phase -- for full Moon.  Sun is at the origin
phi_E = 3.0*math.pi/2.0
phi_M = 3.0*math.pi/2.0

# number of frames
nframes = 500

# number of frames at "freeze" scenes
nfreeze = 100

eps = 1.e-8


def lunar_period():

    # time and timestep
    t = 0
    dt = P_M_syn / (nframes - 1)
    dt = dt*(1 - eps)    # reduce the timestep slightly to make sure
                         # we end at the desired time

    iframe = 0

    # part 0: show the initial configuration
    for n in range(nfreeze):

        draw_earth_moon(t, iframe, connect_EM = 1, label = 1,
                        annotation = "full Moon")
        iframe += 1


    # part I: integrate up to sidereal period
    while (t < P_M_sid):

        draw_earth_moon(t, iframe, connect_EM = 0)

        t += dt
        iframe += 1


    # part II: freeze after 1 lunar sidereal period
    for n in range(nfreeze):

        draw_earth_moon(t, iframe, connect_EM = 0,
                        annotation = "one sidereal lunar orbital period later")

        iframe += 1


    # part III: integrate up to the synodic period -- next full moon
    while (t <= P_M_syn):

        draw_earth_moon(t, iframe, connect_EM = 0)

        t += dt
        iframe += 1


    # part IV: freeze after 1 lunar synodic period
    for n in range(nfreeze):

        draw_earth_moon(t, iframe, connect_EM = 1,
                        annotation = "one synodic period later, full Moon again")

        iframe += 1


def circle(xc, yc, r):
    theta = np.arange(361)*2.0*math.pi/360
    return r*np.cos(theta) + xc, r*np.sin(theta) + yc


def draw_earth_moon(t, iframe, connect_EM = 0, label = 0,
                    annotation = ""):

    plt.clf()

    # compute positions of the Earth
    x_E = x_E_c + r_E*math.cos(omega_E*t + phi_E)
    y_E = y_E_c + r_E*math.sin(omega_E*t + phi_E)


    # compute positions of the Moon
    x_M = x_E + r_M*math.cos(omega_M*t + phi_M)
    y_M = y_E + r_M*math.sin(omega_M*t + phi_M)


    # plot the Sun
    plt.scatter([0], [0], s=1600, marker=(20,1), color="k")
    plt.scatter([0], [0], s=1500, marker=(20,1), color="#FFFF00")

    # plot the Earth
    plt.scatter([x_E], [y_E], s=150, color="b")

    if label:
        plt.text(1.05*x_E + 0.05, 1.05*y_E,
                   "Earth", color="b", fontsize=12)


    # plot the Moon
    plt.scatter([x_M], [y_M], s=50, color="0.5")

    if label:
        plt.text(1.05*x_M + 0.05, 1.05*y_M,
                   "Moon", color="0.5", fontsize=12)


    # plot the orbit of Earth
    xo_E, yo_E = circle(0.0, 0.0, r_E)
    plt.plot(xo_E, yo_E, color="b", linestyle=":", lw=2)


    # plot the orbit of the Moon
    xo_M, yo_M = circle(x_E, y_E, r_M)
    plt.plot(xo_M, yo_M, color="0.5", linestyle=":", lw=2)


    if connect_EM:
        # draw a line connecting the Sun-Earth-Moon to show we are
        # at full Moon again
        plt.plot([0, x_M], [0, y_M], color="0.5", linestyle=":", lw=2)


    if not annotation == "":
        plt.text(-r_M, 1.5*r_M, annotation, fontsize=12)

    plt.axis([-(2*r_M), 1.05*(r_E+r_M),
                -1.05*(r_E+r_M),(2*r_M)])

    plt.axis("off")

    ax = plt.gca()
    ax.set_aspect("equal", "datalim")

    plt.subplots_adjust(left=0.05,right=0.98,bottom=0.05,top=0.98)

    f = plt.gcf()
    f.set_size_inches(7.2,7.2)

    plt.savefig("lunar_period_%04d.png" % iframe)


if __name__== "__main__":
    lunar_period()
