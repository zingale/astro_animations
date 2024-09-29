import matplotlib.pyplot as plt

"""
red black gauss seidel pattern
"""


class Marker:

    def __init__(self, xc, yc, color):

        # a marker is indicated by its center (xc,yc).
        self.xc = xc
        self.yc = yc

        # keep track of the color
        self.color = color


def red_black():

    # define the number of markers in x and y
    nx = 10
    ny = 10

    # define the length of a marker side
    L = 0.8

    # create a list of marker objects, one at each grid location
    markers = []

    for color in [0, 1]:
        for j in range(ny):

            if color == 0:
                ioff = j % 2
            else:
                ioff = 1 - (j % 2)

            for i in range(ioff, nx, 2):
                markers.append(Marker(i, j, color))

    fig, ax = plt.subplots()

    # the margins are funny -- we pick them to ensure that the
    # plot size is an integer multiple of the number of markers in
    # each dimension
    plt.subplots_adjust(left=0.1, right=0.9,
                        bottom=0.1, top=0.9)

    # draw the current state
    for m in markers:

        if (m.color == 1):
            c = "r"
        else:
            c = "k"

        ax.fill([m.xc-L/2, m.xc-L/2,
                 m.xc+L/2, m.xc+L/2,
                 m.xc-L/2],
                [m.yc-L/2, m.yc+L/2,
                 m.yc+L/2, m.yc-L/2,
                 m.yc-L/2],
                color=c)

    ax.axis([-0.5, nx-0.5, -0.5, ny-0.5])
    ax.axis("off")

    fig.set_size_inches(6.0, 6.0)
    fig.tight_layout()

    fig.savefig("rb.png")
    fig.savefig("rb.pdf", bbox_inches="tight", pad_inches=0)


if __name__ == "__main__":
    red_black()
