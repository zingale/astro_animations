import math
import numpy
import pylab
import random

# physical constants
k = 1.3806503e-16    # Boltzmann constant
m_p = 1.67e-24       # proton mass


# a container class for an atom
class atom:

    def __init__ (self, x_init, y_init, vel, angle):
        # angle is assumed to be in radians

        self.x = x_init
        self.y = y_init
        self.u = vel*math.cos(angle)    # x-velocity
        self.v = vel*math.sin(angle)    # y-velocity

    def advance(self, dt, L):
        
        self.x = self.x + self.u*dt
        self.y = self.y + self.v*dt

        if (self.x < 0):
            self.x = 0
            self.u = -self.u

        if (self.x > L):
            self.x = L
            self.u = -self.u

        if (self.y < 0):
            self.y = 0
            self.v = -self.v

        if (self.y > L):
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
    vel = math.sqrt(3.0*k*T/m_p)

    # list to hold the atoms
    atoms = []

    # seed for the random number generator
    seed = 1000
    random.seed(seed)

    # initialize the atoms
    n = 0
    while (n < N):
        x = L*random.random()
        y = L*random.random()
        angle = 2.0*math.pi*random.random()

        atoms.append(atom(x, y, vel, angle))
        
        n += 1


    plot_deltat = 5.e-6
    tmax = 1.e-2

    # compute the timestep 
    dt = 0.1*L/(nstep*vel)


    print "v = ", vel
    print "dt = ", dt
    print "L = ", L

    # take steps, draw the atoms and then move them and write a frame
    iframe = 0
    i = 0
    t = 0
    while (t < tmax):

        pylab.clf()

    
        # draw the box
        pylab.plot([0.0,0.0,L,L,0.0],
                   [0.0,L,L,0.0,0.0],
                   color = "k")


        # draw the atoms
        n = 0
        while (n < len(atoms)):

            pylab.scatter([atoms[n].x], [atoms[n].y])
            n += 1



        # move the atoms
        n = 0
        while (n < len(atoms)):
            
            atoms[n].advance(dt, L)
            n += 1


        # # look for collisions
        # n = 0
        # while (n < len(atoms)):
            
        #     atoms[n].collisionDetect(dt, L, radius)
        #     n += 1


        # axis stuff
        ax = pylab.gca()
        ax.set_aspect("equal", "datalim")

        pylab.axis([0.0-L/20.0,L+L/20.0,
                    0.0-L/20.0,L+L/20.0])
        pylab.subplots_adjust(left=0.05, right=0.95,
                              bottom=0.05, top=0.95)


        pylab.axis("off")

        f = pylab.gcf()
        f.set_size_inches(5.0,5.0)

        if ((t - dt) % plot_deltat > t % plot_deltat):
            print "outputting: ", t, i
            outfile = "thermal_motion_%04d.png" % iframe
            pylab.savefig(outfile)
            iframe += 1
            
        t += dt
        i += 1
        

if __name__== "__main__":
    thermal_motion()

