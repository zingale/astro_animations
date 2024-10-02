import numpy as np
import matplotlib.pyplot as plt
import random

"""simulate radioactive decay by populating a grid with a number of
markers and having each marker have a 50/50 chance of "decaying"
each half-life"""


class Marker:

    def __init__(self, xc, yc):

        # a marker is indicated by its center (xc,yc).
        self.xc = xc
        self.yc = yc

        # the state indicates whether it is the original object (1)
        # or decayed into its daughter product (0)
        self.state = 1

    def decay(self):

        if self.state == 1:

            # random.random() returns a number in the range [0.0, 1.0)
            if random.random() >= 0.5:
                self.state = 0


def num_decayed(markers):
    return len([q for q in markers if q.state == 0])


def radioactive_decay():

    # define the number of markers in x and y
    nx = 50
    ny = 50

    # estimate the number of half-lifes needed to decay all markers
    nest = int(np.log(nx*ny)/np.log(2))

    # allocate storage for the number of markers that are decayed at
    # each half-life.  For safety, we allocate space for 2x our estimate
    hist_decay = np.zeros(2*nest)

    # define the length of a marker side
    L = 0.8

    # create a list of marker objects, one at each grid location
    markers = []
    for i in range(nx):
        for j in range(ny):
            markers.append(Marker(i, j))

    # loop over half-lives, and re-evaluate the marker state
    for t in range(2*nest):

        fig, ax = plt.subplots()

        # the margins are funny -- we pick them to ensure that the
        # plot size is an integer multiple of the number of markers in
        # each dimension
        fig.subplots_adjust(left=0.0493333, right=0.9506666,
                            bottom=0.0493333, top=0.9506666)

        # draw the current state
        for m in markers:

            if m.state == 1:
                c = "r"
            else:
                c = "w"

            ax.fill([m.xc-L/2, m.xc-L/2,
                     m.xc+L/2, m.xc+L/2,
                     m.xc-L/2],
                    [m.yc-L/2, m.yc+L/2,
                     m.yc+L/2, m.yc-L/2,
                     m.yc-L/2],
                    color=c, edgecolor="k")

        nd = num_decayed(markers)
        hist_decay[t] = nd

        ax.axis([-1, nx+1, -1, ny+1])
        ax.axis("off")

        ax.text(-1, -1, f"number decayed = {nd}",
                horizontalalignment="left", verticalalignment="top")

        fig.set_size_inches(7.5, 7.5)

        outfile = f"radioactive_decay_{t:04d}.png"
        fig.savefig(outfile)
        plt.close(fig)

        print(t, nd)

        # if all are decayed, stop making plots
        if nd == nx*ny:
            break

        # now give each marker a chance to "decay"
        for m in markers:
            m.decay()

    # now make a plot of the number that decayed as a function of half-life
    fig, ax = plt.subplots()
    fig.subplots_adjust(left=0.1, right=0.9,
                        bottom=0.1, top=0.9)

    tsmooth = np.arange(1000) * t / 1000.0
    pred = nx*ny*(1.0/2.0)**tsmooth

    ax.plot(tsmooth, pred, color="b")

    time = np.arange(t-1)
    ax.scatter(time, nx*ny-hist_decay[0:t-1], color="b", label="parent")
    ax.scatter(time, hist_decay[0:t-1], color="r", label="daughter")

    ax.axis([0, np.max(time)+1, 0, 1.1*nx*ny])
    ax.set_xlabel("half life")

    leg = ax.legend(loc=2, fontsize="small", frameon=False)

    outfile = f"radioactive_decay_{t+1:04d}.png"
    fig.set_size_inches(7.5, 7.5)
    fig.savefig(outfile)


if __name__ == "__main__":
    radioactive_decay()
