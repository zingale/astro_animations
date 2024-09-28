#!/bin/env python

import random

import numpy as np
import matplotlib.pyplot as plt


class RandomWalk:

    def __init__(self, l, R, seed=1000):

        self.l = l
        self.R = R

        self.seed = seed

        self.max_steps = 5000

        self.x = None
        self.y = None

    def walk(self):

        # set the initial position to be the origin
        x_0 = 0.0
        y_0 = 0.0

        x = []
        y = []

        random.seed(self.seed)

        for _ in range(self.max_steps):

            # compute a random angle
            angle = 2.0*np.pi*random.random()

            # compute the end coordinates of the segment
            x_1 = self.l*np.cos(angle) + x_0
            y_1 = self.l*np.sin(angle) + y_0

            x.append(x_1)
            y.append(y_1)

            # have we hit the edge of our domain?
            if np.sqrt(x_1**2 + y_1**2) >= self.R:
                break

            x_0 = x_1
            y_0 = y_1

        self.x = np.array(x)
        self.y = np.array(y)


def random_walk():

    # length of single step
    l = 1.0

    # maximum radius of domain
    R = 25.0

    # take steps, draw a segment in a random direction, and save the frame
    fig, ax = plt.subplots()
    fig.set_size_inches(7.2, 7.2)

    # draw a circle to indicate the extent of the domain
    npts = 360
    theta = np.arange(npts)*2*np.pi/(npts-1)

    ax.plot(R*np.cos(theta), R*np.sin(theta), color='k')

    plt.subplots_adjust(left=0, right=1.0, bottom=0, top=1.0)

    # do a random walk
    r1 = RandomWalk(l, R, seed=1000)
    r1.walk()

    for n in range(len(r1.x)-1):

        ax.plot([r1.x[n], r1.x[n+1]],
                [r1.y[n], r1.y[n+1]], color='C0')

        ax.axis([-1.1*R, 1.1*R, -1.1*R, 1.1*R])
        ax.axis("off")

        fig.savefig(f"random_walk_{n:04d}.png")

    # now do a second one
    r2 = RandomWalk(l, R, seed=1001)
    r2.walk()

    fig, ax = plt.subplots()
    fig.set_size_inches(7.2, 7.2)

    ax.plot(R*np.cos(theta), R*np.sin(theta), color='k')

    plt.subplots_adjust(left=0, right=1.0, bottom=0, top=1.0)

    # draw the old one lighter
    for n in range(len(r1.x)-1):
        ax.plot([r1.x[n], r1.x[n+1]],
                [r1.y[n], r1.y[n+1]], color='C0', alpha=0.33)

    # now draw the new one
    for n in range(len(r2.x)-1):
        ax.plot([r2.x[n], r2.x[n+1]],
                [r2.y[n], r2.y[n+1]], color='C1')

        ax.axis([-1.1*R, 1.1*R, -1.1*R, 1.1*R])
        ax.axis("off")

        fig.savefig(f"random_walk_{n+len(r1.x):04d}.png")


if __name__== "__main__":
    random_walk()
