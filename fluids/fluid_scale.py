#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from anim_solvers import simple_arrow


class Atom:

    def __init__(self, pos, angle, R=0.5, length=1.0, color="C0"):
        self.pos = pos
        self.angle = angle

        self.R = R
        self.length = length

        self.color = color
        self.circ_theta = np.linspace(0, 2.0*np.pi, 90, endpoint=True)

    def draw(self, ax=None):

        # circle
        xc = self.pos[0] + self.R*np.cos(self.circ_theta)
        yc = self.pos[1] + self.R*np.sin(self.circ_theta)

        ax.fill(xc, yc, color=self.color)

        # arrow
        simple_arrow.draw_arrow(self.pos, self.length, self.angle,
                                color=self.color, ax=ax)


def fluid_element(center, a=1.0, e=0.4, theta=np.pi/4):

    N = 180
    thetas = np.linspace(0, 2.0*np.pi, N, endpoint=True)

    r = a*(1 - e**2)/(1 + e*np.cos(thetas)) + \
        0.025*a*np.cos(7*thetas) + 0.03*a*np.sin(11*thetas)

    x = r*np.cos(thetas)
    y = r*np.sin(thetas)

    # rotate
    xp = center[0] + np.cos(theta)*x - np.sin(theta)*y
    yp = center[1] + np.sin(theta)*x + np.cos(theta)*y
    return xp, yp


def zoom_box(ll, uu, c1s, c1e, c2s, c2e, color="C0", ax=None):

    ax.plot([ll[0], ll[0], uu[0], uu[0], ll[0]],
            [ll[1], uu[1], uu[1], ll[1], ll[1]], color=color)

    ax.plot([c1s[0], c1e[0]], [c1s[1], c1e[1]], color=color)
    ax.plot([c2s[0], c2e[0]], [c2s[1], c2e[1]], color=color)


def draw():

    fig, ax = plt.subplots()

    # box of atoms
    axmin = 0.17
    aymin = 0.2
    axmax = 0.37
    aymax = 0.4

    ax.fill([axmin, axmin, axmax, axmax, axmin],
            [aymin, aymax, aymax, aymin, aymin], color="C8", zorder=-100)

    ax.text(0.5*(axmin+axmax), 0.98*aymin, "atomic scale\n(mean free path)",
            fontsize="large",
            horizontalalignment="center", verticalalignment="top",
            color="black")

    N = 100

    xs = (axmax - axmin) * np.random.random_sample(N) + axmin
    ys = (aymax - aymin) * np.random.random_sample(N) + aymin
    thetas = 2.0*np.pi*np.random.random_sample(N)

    # generate atoms
    atoms = []
    for n in range(100):
        if n % 2 == 0:
            color = "C0"
        else:
            color = "C1"
        atoms.append(Atom((xs[n], ys[n]), thetas[n], R=0.005,
                          length=0.02, color=color))

    for a in atoms:
        a.draw(ax=ax)

    # intermediate scale
    xmin = 0.3
    ymin = 0.5
    xmax = 0.5
    ymax = 0.7

    ax.fill([xmin, xmin, xmax, xmax, xmin],
            [ymin, ymax, ymax, ymin, ymin], color="C8", zorder=-100)

    xe, ye = fluid_element((0.5*(xmin + xmax), 0.5*(ymin + ymax)),
                           a=0.025, e=0.5)
    ax.fill(xe, ye, color="C1")

    ax.text(0.5*(xmin+xmax), 1.02*ymax, "fluid element",
            fontsize="large", horizontalalignment="center", color="black")

    # connect
    zoom_box((0.4, 0.6), (0.404, 0.604),
             (0.4, 0.604), (axmin, aymax),
             (0.404, 0.6), (axmax, aymin), ax=ax)

    # star
    c = (0.8, 0.3)
    R = 0.3

    thetas = np.linspace(0.0, 2.0*np.pi, 180)
    xs = c[0] + R*np.cos(thetas)
    ys = c[1] + R*np.sin(thetas)

    ax.fill(xs, ys, color="C8", zorder=-100)

    # connect
    zoom_box((0.7, 0.2), (0.704, 0.204),
             (0.704, 0.204), (xmax, ymax),
             (0.7, 0.2), (xmin, ymin), color="C2", ax=ax)

    ax.text(1.05, 0.55, "star", horizontalalignment="left",
            fontsize="large", color="black")

    ax.set_aspect("equal", "datalim")
    fig.tight_layout()
    fig.set_size_inches(8, 6.5)

    ax.axis("off")
    fig.savefig("fluid_scale.pdf", dpi=100, bbox_inches="tight", pad_inches=0)
    fig.savefig("fluid_scale.png", dpi=100, transparent=True)


if __name__ == "__main__":
    draw()
