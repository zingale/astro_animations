# a simple class to integrate a satellite around the Earth

import numpy as np
import math

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
        f[2] = -self.GM*y[0]/r**3
        f[3] = -self.GM*y[1]/r**3

        return f
