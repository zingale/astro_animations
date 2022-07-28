import numpy as np
import matplotlib.pyplot as plt


class Wavefront:
    """The wavefront class will hold the location and time that a
    wavefront was emitted

    """

    def __init__(self, x_emit, y_emit, w, t_emit):
        self.x_emit = x_emit
        self.y_emit = y_emit
        self.w = w           # wave propagation speed
        self.t_emit = t_emit


def doppler():

    # emitter velocity (in x-direction)
    vel1 = 1.0
    vel2 = 0.5

    # emitter initial coords
    x_init = 0.0
    y_init = 0.0

    # wave velocity
    w = 2.0

    # wave frequency (# of peaks per second)
    f = 3.0

    # maximum time
    tmax = 10.0
    dt = 0.01

    # create a list of wavefront objects that we can refer to when we
    # want to plot things.  There are f wavefronts emitted per second,
    # so the total number of wavefronts is tmax*f
    t = 0

    wavefronts1 = []
    while t <= tmax:

        x_emit = x_init + vel1 * t
        y_emit = y_init

        wavefronts1.append(Wavefront(x_emit, y_emit, w, t))

        t += 1/f

    t = 0

    wavefronts2 = []
    while t <= tmax:

        x_emit = x_init + vel2*t
        y_emit = y_init

        wavefronts2.append(Wavefront(x_emit, y_emit, w, t))

        t += 1/f

    xmax = x_init + max(vel1, vel2) * tmax

    # we will be drawing circles, so make an array with the polar angle
    npts = 360
    theta = np.linspace(0.0, 2.0 * np.pi, npts, endpoint=True)

    # step forward in time (by dt) and draw any wavefronts that have
    # been emitted
    iframe = 0
    t = 0

    fig, ax = plt.subplots(nrows=2, ncols=1)

    while t <= tmax:

        ax[0].clear()
        ax[1].clear()

        x_source = x_init + vel1*t
        y_source = y_init

        # plot the sources's path
        ax[0].plot([-1.2*xmax, 1.2*xmax], [y_init, y_init], 'k:')

        # draw the source
        ax[0].scatter([x_source], [y_source], color='b')

        # loop over the wavefronts, and draw any that have been
        # emitted so far
        for wf in wavefronts1:

            if wf.t_emit > t:
                break

            r_front = wf.w * (t - wf.t_emit)

            # wavefronts are circles centered on their emitted coordinates
            x_front = wf.x_emit + r_front * np.cos(theta)
            y_front = wf.y_emit + r_front * np.sin(theta)

            ax[0].plot(x_front, y_front, color='r')

        ax[0].plot([-1.2*xmax, 1.2*xmax], [-0.8*xmax, -0.8*xmax], color='k', lw=2)
        ax[0].set_xlim(-1.2*xmax, 1.2*xmax)
        ax[0].set_ylim(-0.8*xmax, 0.8*xmax)
        ax[0].set_axis_off()

        # second source

        x_source = x_init + vel2 * t
        y_source = y_init

        # plot the sources's path
        ax[1].plot([-1.2*xmax, 1.2*xmax], [y_init, y_init], 'k:')

        # draw the source
        ax[1].scatter([x_source], [y_source], color='b')

        # loop over the wavefronts, and draw any that have been
        # emitted so far
        for wf in wavefronts2:

            if wf.t_emit > t:
                break

            r_front = wf.w * (t - wf.t_emit)

            # wavefronts are circles centered on their emitted coordinates
            x_front = wf.x_emit + r_front * np.cos(theta)
            y_front = wf.y_emit + r_front * np.sin(theta)

            ax[1].plot(x_front, y_front, color='g')

        ax[1].plot([-1.2*xmax, 1.2*xmax], [0.8*xmax, 0.8*xmax], color='k', lw=2)
        ax[1].set_xlim(-1.2*xmax, 1.2*xmax)
        ax[1].set_ylim(-0.8*xmax, 0.8*xmax)
        ax[1].set_axis_off()

        fig.subplots_adjust(left=0, right=1.0, bottom=0, top=1.0,
                            wspace=0, hspace=0)

        fig.set_size_inches(4.5, 6.0)

        outfile = f"doppler_{iframe:04d}.png"
        fig.savefig(outfile)

        t += dt
        iframe += 1


if __name__ == "__main__":
    doppler()
