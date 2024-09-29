# integrate a system of binary stars

import numpy as np

# we work in CGS units
G = 6.67428e-8        # cm^3 g^{-1} s^{-2}
M_sun = 1.98892e33    # g
AU = 1.49598e13       # cm
year = 3.1556926e7


class _StarHistory:
    """ a simple container to hold the solution """

    def __init__(self, num_steps):
        self.t = np.zeros(num_steps, np.float64)
        self.x = np.zeros(num_steps, np.float64)
        self.y = np.zeros(num_steps, np.float64)
        self.vx = np.zeros(num_steps, np.float64)
        self.vy = np.zeros(num_steps, np.float64)


# a simple class to integrate a binary system
class Binary:

    def __init__ (self, M1, M2, a, e, theta, annotate=False):
        """
        define a binary system:

        M1 is the mass of object (star/planet) 1
        M2 is the mass of object (star/planet) 2
        a is the sum of the semi-major axes (a1 + a2)
        e is the eccentricity
        theta is an angle to rotate the semi-major axis wrt +x
        """

        self.M1 = M1
        self.M2 = M2
        self.a = a
        self.e = e
        self.theta = theta

        # determine the individual semi-major axes
        # a1 + a2 = a,  M1 a1 = M2 a2
        self.a1 = self.a/(1.0 + self.M1/self.M2)
        self.a2 = self.a - self.a1

        # we put the center of mass at the origin
        # we put star 1 on the -x axis and star 2 on the +x axis
        self.x1_init = -self.a1*(1.0 - self.e)*np.cos(self.theta)
        self.y1_init = -self.a1*(1.0 - self.e)*np.sin(self.theta)

        self.x2_init = self.a2*(1.0 - self.e)*np.cos(self.theta)
        self.y2_init = self.a2*(1.0 - self.e)*np.sin(self.theta)

        # Kepler's laws should tell us the orbital period
        # P^2 = 4 pi^2 (a_star1 + a_star2)^3 / (G (M_star1 + M_star2))
        self.P = np.sqrt(4*np.pi**2*(self.a1 + self.a2)**3/(G*(self.M1 + self.M2)))

        # compute the initial velocities velocities

        # first compute the velocity of the reduced mass at perihelion
        # (C&O Eq. 2.33)
        v_mu = np.sqrt( (G*(self.M1 + self.M2)/(self.a1 + self.a2)) *
                        (1.0 + self.e)/(1.0 - self.e) )

        # then v_star2 = (mu/m_star2)*v_mu
        self.vx2_init = -(self.M1/(self.M1 + self.M2))*v_mu*np.sin(self.theta)
        self.vy2_init =  (self.M1/(self.M1 + self.M2))*v_mu*np.cos(self.theta)

        # then v_star1 = (mu/m_star1)*v_mu
        self.vx1_init =  (self.M2/(self.M1 + self.M2))*v_mu*np.sin(self.theta)
        self.vy1_init = -(self.M2/(self.M1 + self.M2))*v_mu*np.cos(self.theta)

        self.annotate = annotate

        self.orbit1 = None
        self.orbit2 = None


    def integrate(self, dt, tmax):
        """ integrate our system to tmax using a stepsize dt """

        # allocate storage for R-K intermediate results
        # y[0:3] will hold the star1 info, y[4:7] will hold the star2 info
        k1 = np.zeros(8, np.float64)
        k2 = np.zeros(8, np.float64)
        k3 = np.zeros(8, np.float64)
        k4 = np.zeros(8, np.float64)

        y = np.zeros(8, np.float64)

        t = 0.0

        # initial conditions

        # star 1
        y[0] = self.x1_init  # initial x position
        y[1] = self.y1_init  # initial y position

        y[2] = self.vx1_init # initial x-velocity
        y[3] = self.vy1_init # initial y-velocity

        # star 2
        y[4] = self.x2_init  # initial x position
        y[5] = self.y2_init  # initial y position

        y[6] = self.vx2_init # initial x-velocity
        y[7] = self.vy2_init # initial y-velocity


        # how many steps will we need?
        nsteps = int(tmax/dt)
        self.npts = nsteps+1

        # solution storage
        s1 = _StarHistory(self.npts)
        s2 = _StarHistory(self.npts)

        s1.x[0] = self.x1_init
        s1.y[0] = self.y1_init
        s1.vx[0] = self.vx1_init
        s1.vy[0] = self.vy1_init

        s2.x[0] = self.x2_init
        s2.y[0] = self.y2_init
        s2.vx[0] = self.vx2_init
        s2.vy[0] = self.vy2_init

        s1.t[0] = s2.t[0] = t

        for n in range(1, nsteps+1):

            k1[:] = dt*self.rhs(t, y, self.M1, self.M2)
            k2[:] = dt*self.rhs(t+0.5*dt, y[:]+0.5*k1[:], self.M1, self.M2)
            k3[:] = dt*self.rhs(t+0.5*dt, y[:]+0.5*k2[:], self.M1, self.M2)
            k4[:] = dt*self.rhs(t+dt, y[:]+k3[:], self.M1, self.M2)

            y[:] += (1.0/6.0)*(k1[:] + 2.0*k2[:] + 2.0*k3[:] + k4[:])

            t = t + dt

            s1.x[n] = y[0]
            s1.y[n] = y[1]
            s1.vx[n] = y[2]
            s1.vy[n] = y[3]

            s2.x[n] = y[4]
            s2.y[n] = y[5]
            s2.vx[n] = y[6]
            s2.vy[n] = y[7]

            s1.t[n] = s2.t[n] = t

        self.orbit1 = s1
        self.orbit2 = s2

    def kinetic_energies(self):
        KE1 = 0.5 * self.M1 * (self.orbit1.vx**2 + self.orbit1.vy**2)
        KE2 = 0.5 * self.M2 * (self.orbit2.vx**2 + self.orbit2.vy**2)
        return KE1, KE2

    def potential_energy(self):
        PE = -G * self.M1 * self.M2 / np.sqrt((self.orbit1.x - self.orbit2.x)**2 +
                                              (self.orbit1.y - self.orbit2.y)**2)
        return PE

    def rhs(self,t, y, M_star1, M_star2):
        """ the RHS of our system """

        f = np.zeros(8, np.float64)

        # y[0] = x_star1, y[1] = y_star1, y[2] = vx_star1, y[3] = vy_star1
        # y[4] = x_star2, y[5] = y_star2, y[6] = vx_star2, y[7] = vy_star2

        # unpack
        x_star1 = y[0]
        y_star1 = y[1]

        vx_star1 = y[2]
        vy_star1 = y[3]

        x_star2 = y[4]
        y_star2 = y[5]

        vx_star2 = y[6]
        vy_star2 = y[7]


        # distance between stars
        r = np.sqrt((x_star2 - x_star1)**2 + (y_star2 - y_star1)**2)

        f[0] = vx_star1  # d(x_star1) / dt
        f[1] = vy_star1  # d(y_star1) / dt

        f[2] = -G*M_star2*(x_star1 - x_star2)/r**3  # d(vx_star1) / dt
        f[3] = -G*M_star2*(y_star1 - y_star2)/r**3  # d(vy_star1) / dt

        f[4] = vx_star2  # d(x_star2) / dt
        f[5] = vy_star2  # d(y_star2) / dt

        f[6] = -G*M_star1*(x_star2 - x_star1)/r**3  # d(vx_star2) / dt
        f[7] = -G*M_star1*(y_star2 - y_star1)/r**3  # d(vy_star2) / dt

        return f
