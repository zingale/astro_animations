#!/usr/bin/env python

import numpy as np
import matplotlib.pylab as plt

def _rotate(point, center, theta):
    """ apply a rotation matrix

         / cos theta   -sin theta \
         |                        |
         \\ sin theta    cos theta /

        to a point about center, and return the
        transformed point """

    x = point[0] - center[0]
    y = point[1] - center[1]

    return (x*np.cos(theta) - y*np.sin(theta) + center[0],
            x*np.sin(theta) + y*np.cos(theta) + center[1])


def draw_arrow(center, L, rot, frac=0.2, color="k"):
    """ draw an arrow at the position given by center
        of lengtj L and rotated by an angle (radians) rot """

    # our person is 3 segments

    # main line
    main_line_start = (center[0], center[1])
    main_line_end = (center[0] + L, center[1])

    # top part of head
    top_start = (center[0] + L, center[1])
    top_end = (center[0] + (1.0-frac)*L, center[1] + frac*L)

    # bottom part of head
    bottom_start = (center[0] + L, center[1])
    bottom_end = (center[0] + (1.0-frac)*L, center[1] - frac*L)


    # draw our arrow
    lb = _rotate(main_line_start, center, rot)
    le = _rotate(main_line_end, center, rot)

    plt.plot([lb[0], le[0]], [lb[1], le[1]], color=color)

    lb = _rotate(top_start, center, rot)
    le = _rotate(top_end, center, rot)

    plt.plot([lb[0], le[0]], [lb[1], le[1]], color=color)

    lb = _rotate(bottom_start, center, rot)
    le = _rotate(bottom_end, center, rot)

    plt.plot([lb[0], le[0]], [lb[1], le[1]], color=color)


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
        draw_arrow(center, L, np.radians(l - 90), color="r")
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
