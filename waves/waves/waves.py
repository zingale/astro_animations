#!/bin/env python

"""
illustrate the difference between wavelength and frequency
"""

import numpy as np
import matplotlib.pyplot as plt


def waves():

    tmax = 4.0
    nframes = 1200
    dt = tmax/nframes

    lambda_1 = 1.0
    lambda_2 = 0.25

    t = 0.0
    v = 2.0

    xmin = 0
    xmax = 3.0
    npts = 500
    x = np.linspace(xmin, xmax, npts)

    for iframe in range(nframes):

        fig = plt.figure()

        ax = fig.add_subplot(211)

        y = np.sin(2*np.pi*(x - v*t)/lambda_1)

        xpt = 2.0
        ypt = np.sin(2*np.pi*(xpt - v*t)/lambda_1)

        ax.plot(x, y, color="C0")
        ax.scatter([xpt], [ypt],color="C0")

        # draw and annotate a dimension line
        ax.annotate("", (xmin, 1.2), (xmin+lambda_1, 1.2),
                    arrowprops=dict(arrowstyle="|-|", mutation_scale=5.0))

        ax.text(0.5*(xmin + xmin + lambda_1), 1.4,
                rf"$\lambda = $ {lambda_1:3.2f} cm",
                horizontalalignment='center')

        ax.plot([2.0, 2.0], [-1.5, 1.5], "k--")

        ax.axis([xmin, xmax, -1.6, 1.6])
        ax.axis("off")

        ax = fig.add_subplot(212)

        y = np.sin(2*np.pi*(x - v*t)/lambda_2)

        xpt = 2.0
        ypt = np.sin(2*np.pi*(xpt - v*t)/lambda_2)

        ax.plot(x, y, color="C1")
        ax.scatter([xpt], [ypt], color="C1")

        # draw and annotate a dimension line
        ax.annotate("", (xmin, 1.2), (xmin+lambda_2, 1.2),
                    arrowprops=dict(arrowstyle="|-|", mutation_scale=5.0))

        ax.text(0.5*(xmin + xmin + lambda_2), 1.4,
                rf"$\lambda = $ {lambda_2:3.2f} cm",
                horizontalalignment='center')

        ax.plot([2.0, 2.0], [-1.5, 1.5], "k--")

        # draw the velocity vector
        ax.arrow(0.75*xmax, -1.5, 0.13*xmax, 0.0, color="k",
                 length_includes_head=True,
                 head_width=0.15, head_length=0.15, width=0.05)

        ax.text( 0.90*xmax, -1.5,
                 f"v = {v:3.2} cm/s",
                 verticalalignment='center')

        ax.axis([xmin, xmax, -1.6, 1.6])
        ax.axis("off")

        ax.text(0.01, -1.6, f"t = {t:5.3f} s")

        fig.set_size_inches(7.2, 7.2)

        outfile = f"wave_{iframe:04d}.png"
        fig.savefig(outfile)
        plt.close(fig)

        t += dt


if __name__ == "__main__":
    waves()
