#!/bin/env python

import math
import numpy as np
import matplotlib.pylab as plt

# demonstrate synchronous rotation of the Moon

# M. Zingale

# we work in MKS units
G = 6.67428e-11          # [m^3 kg^{-1} s^{-2}]
M_E = 5.9736e24          # [kg]
d = 3.84399e8            # semi-major axis [m]
day = 24.0*60.0*60.0     # [s]


# a simple class to serve as a container for the orbital information
class trajectory:

    def __init__ (self):

        self.npts = -1
        self.maxpoints = 2000
        self.x  = np.zeros(self.maxpoints)
        self.y  = np.zeros(self.maxpoints)
        self.vx = np.zeros(self.maxpoints)
        self.vy = np.zeros(self.maxpoints)
        self.t  = np.zeros(self.maxpoints)

    def integrate(self, x_init, y_init, vx_init, vy_init, dt, tmax):

        # allocate storage for R-K intermediate results
        k1 = np.zeros(4, np.float64)
        k2 = np.zeros(4, np.float64)
        k3 = np.zeros(4, np.float64)
        k4 = np.zeros(4, np.float64)

        y = np.zeros(4, np.float64)
        f = np.zeros(4, np.float64)

        t = 0.0

        # initial conditions
        y[0] = x_init
        y[1] = y_init

        y[2] = vx_init
        y[3] = vy_init

        # store the initial conditions
        self.x[0] = y[0]
        self.y[0] = y[1]

        self.vx[0] = y[2]
        self.vy[0] = y[3]

        self.t[0] = t

        n = 1
        while (n < self.maxpoints and t < tmax):

            f = self.rhs(t, y)
            k1[:] = dt*f[:]

            f = self.rhs(t+0.5*dt, y[:]+0.5*k1[:])
            k2[:] = dt*f[:]

            f = self.rhs(t+0.5*dt, y[:]+0.5*k2[:])
            k3[:] = dt*f[:]

            f = self.rhs(t+dt, y[:]+k3[:])
            k4[:] = dt*f[:]

            y[:] += (1.0/6.0)*(k1[:] + 2.0*k2[:] + 2.0*k3[:] + k4[:])

            t += dt

            self.x[n]  = y[0]
            self.y[n]  = y[1]
            self.vx[n] = y[2]
            self.vy[n] = y[3]
            self.t[n]  = t

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


def doit():

    # set the semi-major axis and eccentricity
    a = 1.0*d
    e = 0.0

    # compute the period of the orbit from Kepler's law and make
    # the timestep by 1/720th of a period
    P_orbital = math.sqrt(4*math.pi*math.pi*a**3/(G*M_E))

    # set the rotation period
    P_rotation = P_orbital

    omega = 2*math.pi/P_rotation

    # set the initial coordinates -- perihelion
    x_init = 0.0
    y_init = -a*(1.0 - e)

    # set the initial velocity
    vx_init = math.sqrt( (G*M_E/a) * (1 + e) / (1 - e))
    vy_init = 0.0


    # the circular orbit will be our baseline
    orbit = trajectory()


    print "period = ", P_orbital/day

    dt = P_orbital/720.0
    tmax = 2*P_orbital

    orbit.integrate(x_init, y_init, vx_init, vy_init, dt, tmax)


    # plotting

    for n in range(orbit.npts):

        plt.clf()

        # plot the Earth
        plt.scatter([0], [0], s=650, color="k")
        plt.scatter([0], [0], s=600, color="b")

        # plot the orbit
        plt.plot(orbit.x[0:orbit.npts], orbit.y[0:orbit.npts], color="0.5", ls=":", lw=2)

        # plot moon (use zorder to put this on top of the orbit line)
        theta = np.arange(180)
        r = 0.05*d  # exaggerate the moon's size
        x_surface = orbit.x[n] + r*np.cos(theta)
        y_surface = orbit.y[n] + r*np.sin(theta)
        plt.fill(x_surface,y_surface,"0.75", edgecolor="0.75", alpha=1.0, zorder=1000)

        # plot a point on the moon's surface
        xpt = orbit.x[n] + r*np.cos(omega*orbit.t[n]+math.pi/2.0)
        ypt = orbit.y[n] + r*np.sin(omega*orbit.t[n]+math.pi/2.0)
        plt.scatter([xpt],[ypt],s=25,color="k")


        plt.axis([-1.1*d,1.1*d,-1.1*d,1.1*d])
        plt.axis("off")

        plt.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

        ax = plt.gca()
        ax.set_aspect("equal", "datalim")


        f = plt.gcf()
        f.set_size_inches(7.2,7.2)

        plt.savefig("moon_rotation_%04d.png" % n)


if __name__== "__main__":
    doit()
