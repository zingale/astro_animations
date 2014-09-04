import math
import numpy as np
import matplotlib.pyplot as plt

# Shoot a projective horizontally some distance above Earth at various
# speeds and watch the resulting orbit
#
# In this version, we just consider orbits with the horizontal 
# velocity <= the circular velocity.


# M. Zingale (2008-09-14)

# we work in CGS units
G = 6.67e-8
M_E = 5.9742e27
R_E = 6.378e8


# a simple class to serve as a container for the orbital information
class trajectory:
    
    def __init__ (self, GM=4*np.pi**2, R_crash=1.0):

        self.npts = -1

        self.maxpoints = 2000
        self.x = np.zeros(self.maxpoints)
        self.y = np.zeros(self.maxpoints)
        self.t = np.zeros(self.maxpoints)

        self.GM = GM
        self.R_crash = R_crash


    def integrate(self, vinit, hinit, dt, max_rad):

        SMALL = 1.e-16

        # allocate storage for R-K intermediate results
        k1 = np.zeros(4, np.float64)
        k2 = np.zeros(4, np.float64)
        k3 = np.zeros(4, np.float64)
        k4 = np.zeros(4, np.float64)


        y = np.zeros(4, np.float64)
        f = np.zeros(4, np.float64)

        t = 0.0

        # initial conditions
        y[0] = 0       # initial x position
        y[1] = hinit   # initial height

        y[2] = vinit   # initial tangental velocity
        y[3] = 0.0     # initial vertical (radial) velocity

        self.x[0] = y[0]
        self.y[0] = y[1]
        self.t[0] = t

        r = np.sqrt(y[0]**2 + y[1]**2)

        numorbits = 0
        angle_old = np.pi/2   # we are initially vertical

        n = 1
        while (n < self.maxpoints and r > self.R_crash and r <= max_rad and 
               numorbits == 0):

            k1[:] = dt*self.rhs(t, y)
            k2[:] = dt*self.rhs(t+0.5*dt, y[:]+0.5*k1[:])
            k3[:] = dt*self.rhs(t+0.5*dt, y[:]+0.5*k2[:])
            k4[:] = dt*self.rhs(t+dt, y[:]+k3[:])

            y[:] += (1.0/6.0)*(k1[:] + 2.0*k2[:] + 2.0*k3[:] + k4[:])

            t = t + dt

            self.x[n] = y[0]
            self.y[n] = y[1]
            self.t[n] = t

            # compute the radius to decide if we collided with Earth
            r = np.sqrt(y[0]**2 + y[1]**2)

            # compute the angle wrt to vertical and compare to the 
            # previous one.  If we go from > pi/2 to <= pi/2 then we have
            # completed an orbit
            angle = math.atan2(y[1],(y[0] + SMALL))

            if angle_old > np.pi/2 and angle <= np.pi/2:
                numorbits = 1

            angle_old = angle

            n += 1

    
        self.npts = n


    def rhs(self, t, y):

        f = np.zeros(4, np.float64)

        # y[0] = x, y[1] = y, y[2] = v_x, y[3] = v_y
        f[0] = y[2]
        f[1] = y[3]
        r = np.sqrt(y[0]**2 + y[1]**2)
        f[2] = -G*M_E*y[0]/r**3
        f[3] = -G*M_E*y[1]/r**3

        return f
    


def orbit():

    # set the height from which we launch things
    hinit = 1.5*R_E
    max_rad = 25*R_E

    # the circular orbit will be our baseline
    orbit_circ = trajectory(GM=G*M_E, R_crash=R_E)

    vinit_circ = np.sqrt(G*M_E/hinit)
    tmax_circ = 2*np.pi*hinit/vinit_circ

    # we want dt to be constant for all of our orbits so we can compare
    # the velocities in the animation
    dt = tmax_circ/720.0        
    
    orbit_circ.integrate(vinit_circ, hinit, dt, max_rad)
    print "circ: ", orbit_circ.npts

    # orbit 1 (0.25 v_c)
    orbit_1 = trajectory(GM=G*M_E, R_crash=R_E)
    vinit_1 = 0.25*vinit_circ

    orbit_1.integrate(vinit_1, hinit, dt, max_rad)
    print "1: ", orbit_1.npts

    # orbit 2 (0.5 v_c)
    orbit_2 = trajectory(GM=G*M_E, R_crash=R_E)
    vinit_2 = 0.5*vinit_circ

    orbit_2.integrate(vinit_2, hinit, dt, max_rad)
    print "2: ", orbit_2.npts

    # orbit 3 (0.75 v_c)
    orbit_3 = trajectory(GM=G*M_E, R_crash=R_E)
    vinit_3 = 0.75*vinit_circ

    orbit_3.integrate(vinit_3, hinit, dt, max_rad)
    print "3: ", orbit_3.npts



    # ================================================================
    # plotting
    # ================================================================

    img = plt.imread("earth.png")

    # plot the orbits one by one
    iframe = 0

    # v1
    for n in range(orbit_1.npts):

        plt.clf()
        plt.title(r"$v = 0.25 v_\mathrm{circular}$", color="g", fontsize=20)

        # draw Earth -- we will use units that are in terms of Earth radii
        plt.imshow(img,extent = [-R_E,R_E,-R_E,R_E])

        # plot the current orbit
        plt.plot(orbit_1.x[0:n],orbit_1.y[0:n], color="g")


        plt.axis([-3*R_E,3*R_E,-3*R_E,2*R_E])
        
        plt.axis("off")
        
        ax = plt.gca()
        ax.set_aspect("equal", "datalim")

        f = plt.gcf()
        f.set_size_inches(9.6,7.2)


        plt.text(0.05, 0.06, "Earth image credit:", transform=f.transFigure,
                 fontsize=7, color="0.50")
        plt.text(0.05, 0.04, "NASA/Apollo 17", transform=f.transFigure,
                    fontsize=7, color="0.50")

        # print the time
        plt.text(0.7,0.1, "time = %6.4f hrs." % (orbit_1.t[n]/3600.),
                 transform=f.transFigure)

        plt.savefig("achieveorbit_%04d.png" % iframe)

        iframe += 1


    # v2
    for n in range(orbit_2.npts):

        plt.clf()
        plt.title(r"$v = 0.5 v_\mathrm{circular}$", color="r", fontsize=20)

        # draw Earth
        plt.imshow(img,extent = [-R_E,R_E,-R_E,R_E])


        # plot the previous orbit
        plt.plot(orbit_1.x[0:orbit_1.npts-1],
                   orbit_1.y[0:orbit_1.npts-1], color="g", alpha=0.33)

        # plot the current orbit
        plt.plot(orbit_2.x[0:n],orbit_2.y[0:n], color="r")

        plt.axis([-3*R_E,3*R_E,-3*R_E,2*R_E])
        
        plt.axis("off")

        ax = plt.gca()
        ax.set_aspect("equal", "datalim")

        f = plt.gcf()
        f.set_size_inches(9.6,7.2)

        plt.text(0.05, 0.06, "Earth image credit:", transform=f.transFigure,
                 fontsize=7, color="0.50")
        plt.text(0.05, 0.04, "NASA/Apollo 17", transform=f.transFigure,
                    fontsize=7, color="0.50")

        # print the time
        plt.text(0.7,0.1, "time = %6.4f hrs." % (orbit_2.t[n]/3600.),
                 transform=f.transFigure)

        plt.savefig("achieveorbit_%04d.png" % iframe)

        iframe += 1


    # v3
    for n in range(orbit_3.npts):

        plt.clf()
        plt.title(r"$v = 0.75 v_\mathrm{circular}$", 
                    color="b", fontsize=20)

        # draw Earth
        plt.imshow(img,extent = [-R_E,R_E,-R_E,R_E])

        # plot the previous orbits
        plt.plot(orbit_1.x[0:orbit_1.npts-1],
                   orbit_1.y[0:orbit_1.npts-1], color="g", alpha=0.33)
        plt.plot(orbit_2.x[0:orbit_2.npts-1],
                   orbit_2.y[0:orbit_2.npts-1], color="r", alpha=0.33)

        # plot the current orbit
        plt.plot(orbit_3.x[0:n],orbit_3.y[0:n], color="b")

        plt.axis([-3*R_E,3*R_E,-3*R_E,2*R_E])
        
        plt.axis("off")

        ax = plt.gca()
        ax.set_aspect("equal", "datalim")

        f = plt.gcf()
        f.set_size_inches(9.6,7.2)

        plt.text(0.05, 0.06, "Earth image credit:", transform=f.transFigure,
                 fontsize=7, color="0.50")
        plt.text(0.05, 0.04, "NASA/Apollo 17", transform=f.transFigure,
                    fontsize=7, color="0.50")

        # print the time
        plt.text(0.7,0.1, "time = %6.4f hrs." % (orbit_3.t[n]/3600.),
                 transform=f.transFigure)

        plt.savefig("achieveorbit_%04d.png" % iframe)

        iframe += 1


    # v = v_c
    for n in range(orbit_circ.npts):

        plt.clf()
        plt.title(r"$v = v_\mathrm{circular}$", color="k", fontsize=20)

        # draw Earth
        plt.imshow(img,extent = [-R_E,R_E,-R_E,R_E])

        # plot the previous orbits
        plt.plot(orbit_1.x[0:orbit_1.npts-1],
                   orbit_1.y[0:orbit_1.npts-1], color="g", alpha=0.33)
        plt.plot(orbit_2.x[0:orbit_2.npts-1],
                   orbit_2.y[0:orbit_2.npts-1], color="r", alpha=0.33)
        plt.plot(orbit_3.x[0:orbit_3.npts-1],
                   orbit_3.y[0:orbit_3.npts-1], color="b", alpha=0.33)
        
        # plot the current orbit
        plt.plot(orbit_circ.x[0:n+1],orbit_circ.y[0:n+1],color="k")

        plt.axis([-3*R_E,3*R_E,-3*R_E,2*R_E])
        
        plt.axis("off")

        ax = plt.gca()
        ax.set_aspect("equal", "datalim")

        f = plt.gcf()
        f.set_size_inches(9.6,7.2)

        plt.text(0.05, 0.06, "Earth image credit:", transform=f.transFigure,
                 fontsize=7, color="0.50")
        plt.text(0.05, 0.04, "NASA/Apollo 17", transform=f.transFigure,
                    fontsize=7, color="0.50")

        # print the time
        plt.text(0.7,0.1, "time = %6.4f hrs." % (orbit_circ.t[n]/3600.),
                 transform=f.transFigure)

        plt.savefig("achieveorbit_%04d.png" % iframe)

        iframe += 1

    
if __name__== "__main__":
    orbit()


    
        
