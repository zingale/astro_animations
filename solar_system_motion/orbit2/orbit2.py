import math
import numpy
import pylab

# compute the orbit of a 2 planets around the Sun using fourth-order accurate
# Runge-Kutta

# M. Zingale (2008-09-01)

# we work in units of AU, yr, and M_sun
# in these units, G = 4 pi^2

# physical constants
pi = math.pi
G = 4.0*pi*pi  # AU^3 M_sun^{-1} yr^{-2}
M_sun = 1.0

# planet data
ecc_A = 0.0   # eccentricity of planet A
ecc_B = 0.4   # eccecntricity of planet B

a_A = 1.0     # semi-major axis of planet A
a_B = 1.5874  # semi-major axis of planet B

# integration data
nsteps_year = 365   # number of steps per year
nyears = 4          # total integration time (years)


def integrate():

    # allocate storage for R-K intermediate results
    k1 = numpy.zeros(8, numpy.float64)
    k2 = numpy.zeros(8, numpy.float64)
    k3 = numpy.zeros(8, numpy.float64)
    k4 = numpy.zeros(8, numpy.float64)

    # y[0], y[1] are (x,y); y[2], y[3] are (v_x, v_y) for planet 1
    # y[4], y[5] are (x,y); y[6], y[7] are (v_x, v_y) for planet 2
    y = numpy.zeros(8, numpy.float64)
    f = numpy.zeros(8, numpy.float64)

    t = 0.0
    dt = 0.0


    # ================================================================
    # planet A initialization
    # ================================================================

    # our coordinate system puts the primary foci at the origin.

    # set the initial conditions.  The initial position is perihelion
    y[0] = a_A*(1.0 - ecc_A)   # x 
    y[1] = 0.0                 # y

    # at perihelion, all the veloicity is in the y-direction.
    y[2] = 0.0      # v_x

    # v_y^2 = (GM/a) (1+e)/(1-e)  (see C&O Eq. 2.33 for example)
    # This is the perihelion velocity.  For a = 1, e = 0, v = 2 pi
    y[3] = math.sqrt( (G*M_sun/a_A) * (1.0 + ecc_A) / (1.0 - ecc_A))


    # ================================================================
    # planet B initialization
    # ================================================================

    # our coordinate system puts the primary foci at the origin.

    # set the initial conditions.  The initial position is aphelion
    y[4] = -a_B*(1.0 + ecc_B)  # x 
    y[5] = 0.0                 # y

    # at aphelion, all the veloicity is in the negative y-direction.
    y[6] = 0.0      # v_x

    # v_y^2 = (GM/a) (1-e)/(1+e)  (see C&O Eq. 2.33 for example)
    # This is the aphelion velocity.  For a = 1, e = 0, v = 2 pi
    y[7] = -math.sqrt( (G*M_sun/a_B) * (1.0 - ecc_B) / (1.0 + ecc_B))



    # set the initial timestep
    dt = 1.0/nsteps_year

    # compute the total number of steps needed 
    nsteps = int(nyears*nsteps_year)


    # allocate storage to save the positions for all steps
    x_pos_A = numpy.zeros(nsteps, numpy.float64)
    y_pos_A = numpy.zeros(nsteps, numpy.float64)

    x_pos_B = numpy.zeros(nsteps, numpy.float64)
    y_pos_B = numpy.zeros(nsteps, numpy.float64)

    time  = numpy.zeros(nsteps, numpy.float64)


    n = 0  
    while (n < nsteps):

        f = rhs(t, y)
        k1[:] = dt*f[:]

        f = rhs(t+0.5*dt, y[:]+0.5*k1[:])
        k2[:] = dt*f[:]

        f = rhs(t+0.5*dt, y[:]+0.5*k2[:])
        k3[:] = dt*f[:]

        f = rhs(t+dt, y[:]+k3[:])
        k4[:] = dt*f[:]

        y[:] += (1.0/6.0)*(k1[:] + 2.0*k2[:] + 2.0*k3[:] + k4[:])

        # save the data for plotting later
        x_pos_A[n] = y[0]
        y_pos_A[n] = y[1]

        x_pos_B[n] = y[4]
        y_pos_B[n] = y[5]

        time[n]  = t

        t = t + dt
        n += 1


    # ================================================================
    # plotting
    # ================================================================

    # turn on interactive mode 
    # pylab.ion()

    n = 0
    while (n < nsteps):
        pylab.clf()

        # plot the foci
        pylab.scatter([0],[0],s=250,marker=(5,1),color="k")
        pylab.scatter([0],[0],s=200,marker=(5,1),color="y")

        # plot planet A
        pylab.plot(x_pos_A,y_pos_A,color="r")
        pylab.scatter([x_pos_A[n]],[y_pos_A[n]],s=100,color="r")

        # plot planet B
        pylab.plot(x_pos_B,y_pos_B,color="b")
        pylab.scatter([x_pos_B[n]],[y_pos_B[n]],s=100,color="b")

        pylab.axis([-2.5,1.5,-2,2])

        ax = pylab.gca()
        ax.set_aspect("equal", "datalim")

        f = pylab.gcf()
        f.set_size_inches(6.0,6.0)

        pylab.xlabel("AU")
        pylab.ylabel("AU")
        pylab.text(-2.0,-1.75, "time = %6.3f yr" % time[n])
        pylab.text(-2.0,1.75,"a = %6.3f, e = %5.2f" % (a_A, ecc_A) ,color="r")
        pylab.text(-2.0,1.55, "a = %6.3f, e = %5.2f" % (a_B, ecc_B) ,color="b")
        pylab.title("Planetary Orbits")

        outfile = "frame_orbit_%04d.png" % n
        pylab.savefig(outfile)
        n += 1



def rhs(t,y):

    f = numpy.zeros(8, numpy.float64)

    # ================================================================
    # planet A RHS
    # ================================================================

    # y[0] = x, y[1] = y, y[2] = v_x, y[3] = v_y
    f[0] = y[2]
    f[1] = y[3]
    r = numpy.sqrt(y[0]**2 + y[1]**2)
    f[2] = -G*M_sun*y[0]/r**3
    f[3] = -G*M_sun*y[1]/r**3


    # ================================================================
    # planet B RHS
    # ================================================================

    # y[4] = x, y[5] = y, y[6] = v_x, y[7] = v_y
    f[4] = y[6]
    f[5] = y[7]
    r = numpy.sqrt(y[4]**2 + y[5]**2)
    f[6] = -G*M_sun*y[4]/r**3
    f[7] = -G*M_sun*y[5]/r**3

    return f


    
if __name__== "__main__":
    integrate()


    
        
