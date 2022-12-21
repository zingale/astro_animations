# reconstruct - evolve - average: demonstrate what happens when we don't
# limit

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import grid_plot as gp

# font sizes
mpl.rcParams['font.size'] = 20
mpl.rcParams['legend.fontsize'] = 'large'
mpl.rcParams['figure.titlesize'] = 'large'

def evolve(gr, pl, C, num, nolimit=1):

    #-------------------------------------------------------------------------
    # first frame -- the original cell-averages

    nzones = gr.nx

    plt.clf()

    gr.draw_grid()

    labels = ["$i-2$", "$i-1$", "$i$", "$i+1$", "$i+2$"]
    indices = [gr.ng+nzones//2-2, gr.ng+nzones//2-1,
               gr.ng+nzones//2, gr.ng+nzones//2+1, gr.ng+nzones//2+2]

    for i, l in zip(indices, labels):
        gr.label_center(i, l)

    # draw cell averages
    for n in range(gr.ilo, gr.ihi+1):
        pl.draw_cell_avg(n, color="r")

    gr.clean_axes(ylim=(-0.75, 2.5))

    ax = plt.gca()

    plt.text(0.5, 0.85, "initial state (cell averages)",
             horizontalalignment="center",
             fontsize="large", color="b", transform=ax.transAxes)

    plt.text(0.5, 0.95, "Piecewise Linear Method for Linear Advection",
             horizontalalignment="center",
             fontsize="x-large", color="k", transform=ax.transAxes)

    f = plt.gcf()
    f.set_size_inches(12.8,7.2)

    if (nolimit):
        plt.savefig("rea-nolimit-start_%3.3d.png" % (num), dpi=150)
    else:
        plt.savefig("rea-start_%3.3d.png" % (num), dpi=150)

    #-------------------------------------------------------------------------
    # second frame -- reconstruction

    # draw
    plt.clf()

    gr.draw_grid()

    for i, l in zip(indices, labels):
        gr.label_center(i, l)

    # draw cell averages and slopes
    for n in range(gr.ilo, gr.ihi+1):
        pl.draw_cell_avg(n, color="0.5", ls=":")
        pl.draw_slope(n, color="r")

    gr.clean_axes(ylim=(-0.75, 2.5))

    plt.text(0.5, 0.85, "reconstructed slopes",
             horizontalalignment="center",
             fontsize="large", color="b", transform=ax.transAxes)

    plt.text(0.5, 0.95, "Piecewise Linear Method for Linear Advection",
             horizontalalignment="center",
             fontsize="x-large", color="k", transform=ax.transAxes)

    f = plt.gcf()
    f.set_size_inches(12.8,7.2)

    if (nolimit):
        plt.savefig("rea-nolimit-reconstruction_%3.3d.png" % (num), dpi=150)
    else:
        plt.savefig("rea-reconstruction_%3.3d.png" % (num), dpi=150)


    #-------------------------------------------------------------------------
    # third frame -- evolve

    # draw

    plt.clf()

    gr.draw_grid()

    for i, l in zip(indices, labels):
        gr.label_center(i, l)

    # draw cell averages and slopes
    for n in range(gr.ilo, gr.ihi+1):
        pl.draw_slope(n, color="0.5", ls=":")

    # evolve
    for n in range(gr.ilo, gr.ihi+1):
        pl.evolve_to_right(n, C, color="r")

    gr.clean_axes(ylim=(-0.75, 2.5))

    plt.text(0.5, 0.85, f"evolved with C = {C}",
             horizontalalignment="center",
             fontsize="large", color="b", transform=ax.transAxes)

    plt.text(0.5, 0.95, "Piecewise Linear Method for Linear Advection",
             horizontalalignment="center",
             fontsize="x-large", color="k", transform=ax.transAxes)

    f = plt.gcf()
    f.set_size_inches(12.8,7.2)

    if (nolimit):
        plt.savefig("rea-nolimit-evolve_%3.3d.png" % (num), dpi=150)
    else:
        plt.savefig("rea-evolve_%3.3d.png" % (num), dpi=150)


    #-------------------------------------------------------------------------
    # fourth frame -- re-average

    # left states (we don't need the right state when u > 0)
    al = gr.scratch_array()

    for n in range(gr.ilo, gr.ihi+2):
        al[n] = pl.a[n-1] + 0.5*(1 - C)*pl.slope[n-1]

    # the Riemann problem just picks the right state.  Do a conservative
    # update
    anew = gr.scratch_array()

    anew[gr.ilo:gr.ihi+1] = pl.a[gr.ilo:gr.ihi+1] + \
        C*(al[gr.ilo:gr.ihi+1] - al[gr.ilo+1:gr.ihi+2])

    plt.clf()

    gr.draw_grid()

    for i, l in zip(indices, labels):
        gr.label_center(i, l)

    # show the evolved profiles from the old time
    for n in range(gr.ilo, gr.ihi+1):
        pl.evolve_to_right(n, C, color="0.5", ls=":")

    pl.a[:] = anew[:]
    pl.fill_zero_gradient()
    pl.calculate_slopes()

    for n in range(pl.gr.ilo, pl.gr.ihi+1):
        pl.draw_cell_avg(n, color="r")

    gr.clean_axes(ylim=(-0.75, 2.5))

    plt.text(0.5, 0.85, "averaged profile (final state)",
             horizontalalignment="center",
             fontsize="large", color="b", transform=ax.transAxes)

    plt.text(0.5, 0.95, "Piecewise Linear Method for Linear Advection",
             horizontalalignment="center",
             fontsize="x-large", color="k", transform=ax.transAxes)

    f = plt.gcf()
    f.set_size_inches(12.8,7.2)

    if (nolimit):
        plt.savefig("rea-nolimit-final_%3.3d.png" % (num), dpi=150)
    else:
        plt.savefig("rea-final_%3.3d.png" % (num), dpi=150)

    return anew


def main():

    ainit = np.array([1.0, 1.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
    nzones = len(ainit)

    nolimit = 1

    # CFL number
    C = 0.7

    gr = gp.FVGrid(nzones, ng=4)

    a = gr.scratch_array()
    a[gr.ilo:gr.ihi+1] = ainit[:]

    pl = gp.PiecewiseLinear(gr, a, nolimit=nolimit)

    # loop
    for i in range(1,9):
        pl.fill_zero_gradient()
        evolve(gr, pl, C, i, nolimit=nolimit)



if __name__ == "__main__":
    main()
