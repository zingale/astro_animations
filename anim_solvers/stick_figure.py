#!/usr/bin/env python

import numpy as np
import matplotlib.pylab as plt

def _rotate(point, center, theta):
    """ apply a rotation matrix

         / cos theta   -sin theta \
         |                        |
         \ sin theta    cos theta /

        to a point about center, and return the
        transformed point """

    x = point[0] - center[0]
    y = point[1] - center[1]

    return (x*np.cos(theta) - y*np.sin(theta) + center[0],
            x*np.sin(theta) + y*np.cos(theta) + center[1])


def draw_person(center, L, rot, color="k"):
    """ draw a stick figure at the position given by center
        of height L and rotated by an angle (radians) rot """

    # our person is 4 segments

    # head: upper L/4, so the center of the head is L/8
    # above bottom
    head_center = (center[0], center[1]+0.375*L)
    head_radius = L/8.0

    # arms start at the center
    left_arm_end = (center[0]-0.25*L, center[1]+0.25*L)
    right_arm_end = (center[0]+0.25*L, center[1]+0.25*L)

    # torso is the middle two segments
    torso_start = (center[0], center[1]-0.25*L)
    torso_end = (center[0], center[1]+0.25*L)

    # feet start at torso_star
    left_foot_end = (center[0]-0.2*L, center[1]-0.5*L)
    right_foot_end = (center[0]+0.2*L, center[1]-0.5*L)

    # draw our person
    theta = np.radians(np.arange(0,361))
    hc = _rotate(head_center, center, rot)

    plt.fill(hc[0] + head_radius*np.cos(theta),
             hc[1] + head_radius*np.sin(theta), color=color)

    cc = _rotate(center, center, rot)

    lc = _rotate(left_arm_end, center, rot)
    rc = _rotate(right_arm_end, center, rot)

    plt.plot([cc[0], lc[0]], [cc[1], lc[1]], color=color)
    plt.plot([cc[0], rc[0]], [cc[1], rc[1]], color=color)


    ts = _rotate(torso_start, center, rot)
    te = _rotate(torso_end, center, rot)

    plt.plot([ts[0], te[0]], [ts[1], te[1]], color=color)


    lf = _rotate(left_foot_end, center, rot)
    rf = _rotate(right_foot_end, center, rot)

    plt.plot([ts[0], lf[0]], [ts[1], lf[1]], color=color)
    plt.plot([ts[0], rf[0]], [ts[1], rf[1]], color=color)


def doit():
    # test it out

    L = 1

    R = 10

    theta = np.radians(np.arange(0,361))

    # draw a circle
    plt.plot(R*np.cos(theta), R*np.sin(theta), c="b")


    # draw some people
    angles = [30, 60, 90, 120, 180, 270, 300]

    for l in angles:
        center = ( (R + 0.5*L)*np.cos(np.radians(l)),
                   (R + 0.5*L)*np.sin(np.radians(l)) )
        draw_person(center, L, np.radians(l - 90), color="r")
        L = 1.1*L

    plt.axis("off")

    ax = plt.gca()
    ax.set_aspect("equal", "datalim")


    plt.subplots_adjust(left=0.05, right=0.98, bottom=0.05, top=0.98)
    plt.axis([-1.2*R, 1.2*R, -1.2*R, 1.2*R])

    f = plt.gcf()
    f.set_size_inches(6.0, 6.0)

    plt.savefig("test.png")


if __name__ == "__main__":
    doit()
