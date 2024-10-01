import numpy as np
import matplotlib.pyplot as plt
import random

# physical constants
k = 1.3806503e-16    # Boltzmann constant
m_p = 1.67e-24       # proton mass


# a container class for an atom
class Atom:

    def __init__(self, x_init, y_init, vel, angle):
        # angle is assumed to be in radians

        self.x = x_init
        self.y = y_init
        self.u = vel*np.cos(angle)    # x-velocity
        self.v = vel*np.sin(angle)    # y-velocity

    def advance(self, dt, L):

        self.x = self.x + self.u*dt
        self.y = self.y + self.v*dt

        if self.x < 0:
            self.x = 0
            self.u = -self.u
        elif self.x > L:
            self.x = L
            self.u = -self.u

        if self.y < 0:
            self.y = 0
            self.v = -self.v
        elif self.y > L:
            self.y = L
            self.v = -self.v


def thermal_motion():

    # length of box size (cm)
    L = 100.0

    # minimum number of timesteps needed to cross the box
    nstep = 100

    # number of box crossings to model
    ncross = 10

    # number of atoms
    N = 20

    # temperature
    #T = 100.0
    T = 1000.0
    vel = np.sqrt(3.0*k*T/m_p)

    # list to hold the atoms
    atoms = []

    # seed for the random number generator
    seed = 1000
    random.seed(seed)

    # initialize the atoms
    for n in range(N):
        x = L * random.random()
        y = L * random.random()
        angle = 2.0*np.pi*random.random()

        atoms.append(Atom(x, y, vel, angle))

    plot_deltat = 5.e-6
    tmax = 1.e-2

    # compute the timestep
    dt = 0.1*L/(nstep*vel)

    # take steps, draw the atoms and then move them and write a frame
    iframe = 0
    t = 0
    while t < tmax:

        do_plot = (t - dt) % plot_deltat > t % plot_deltat

        if do_plot:
            fig, ax = plt.subplots()

            # draw the box
            ax.plot([0.0, 0.0, L, L, 0.0],
                    [0.0, L, L, 0.0, 0.0],
                    color="k")

            # draw the atoms
            for n in range(len(atoms)):
                ax.scatter([atoms[n].x], [atoms[n].y], color="C0")

            ax.set_aspect("equal", "datalim")
            ax.axis([0.0-L/20.0, L+L/20.0,
                     0.0-L/20.0, L+L/20.0])
            ax.axis("off")

            fig.subplots_adjust(left=0.05, right=0.95,
                                bottom=0.05, top=0.95)

            fig.set_size_inches(7.2, 7.2)

            outfile = f"thermal_motion_{iframe:04d}.png"
            fig.savefig(outfile)
            plt.close(fig)

            iframe += 1

        # move the atoms
        for n in range(len(atoms)):
            atoms[n].advance(dt, L)

        t += dt


if __name__ == "__main__":
    thermal_motion()
