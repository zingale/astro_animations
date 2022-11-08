import numpy as np
import matplotlib.pyplot as plt


def coulomb_barrier():

    fig, ax = plt.subplots()

    npts = 200
    r = np.arange(npts, dtype=np.float64) / npts

    SMALL = 1.e-16

    # Coulomb potential
    C = 1/(r + SMALL)

    # interior, strong potential
    V = np.zeros(npts) - 2.0

    # merge the two
    U = np.where(C < 10.0, C, V)

    ax.plot(r, U, 'r', linewidth=2)

    # draw the axes manually
    ax.axis("off")

    # y-axis
    ax.arrow(0, 1.1*min(U), 0, 1.1*max(U) - 1.1*min(U),
             fc='k', head_width=0.015, head_length=0.18)

    ax.text(-0.05*max(r), 1.1*max(U), r"$U$", color="k", size=15)

    # x-axis
    ax.arrow(0, 0, 1.05*max(r), 0,
             fc='k', head_width=0.12, head_length=0.02)

    ax.text(1.05*max(r), -0.025*max(U), r"$r$", color="k", size=15,
            verticalalignment="top")

    # draw a line representing a classical particle's energy
    E = 0.3*max(U)
    ax.plot([0, max(r)], [E, E], 'b--')
    ax.text(0.75*max(r), 1.05*E, "incoming proton KE", color='b',
            horizontalalignment="center")

    # find the x-coord of the classical turning point
    x = 1.0 / E
    arrow_params = {'length_includes_head': True}

    ax.arrow(x, 0, 0, E, fc='0.5', ec='0.5',
             head_width=0.015, head_length=0.25,
             **arrow_params)

    ax.text(x, -0.025*max(U), r"$r_0$", color="0.5", size=15,
            horizontalalignment="center", verticalalignment="top")

    fig.set_size_inches(6.0, 6.0)
    fig.tight_layout()
    fig.savefig("coulomb_barrier.png", dpi=120)


if __name__ == "__main__":
    coulomb_barrier()
