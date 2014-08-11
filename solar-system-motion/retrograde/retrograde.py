import math
import numpy
import pylab

# compute the orbit of Earth and Mars around the Sun using fourth-order 
# accurate Runge-Kutta.  Include a line connecting the two to show the
# line of sight, and therefore retrograde motion.

# M. Zingale (2008-09-02)

# we work in units of AU, yr, and M_sun
# in these units, G = 4 pi^2

# physical constants
pi = math.pi
G = 4.0*pi*pi  # AU^3 M_sun^{-1} yr^{-2}
M_sun = 1.0

# planet data
ecc_E = 0.016710  # eccentricity of planet Earth
ecc_M = 0.093315  # eccecntricity of planet Mars

a_E = 1.0       # semi-major axis of planet Earth
a_M = 1.523679  # semi-major axis of planet Mars

# integration data
nsteps_year = 365   # number of steps per year
nyears = 0.6        # total integration time (years)


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
    # Earth initialization
    # ================================================================

    # our coordinate system puts the primary foci at the origin.

    # set the initial conditions.  The initial position is perihelion
    #y[0] = a_E*(1.0 - ecc_E)   # x 
    #y[1] = 0.0                 # y

    # at perihelion, all the veloicity is in the y-direction.
    #y[2] = 0.0      # v_x

    # v_y^2 = (GM/a) (1+e)/(1-e)  (see C&O Eq. 2.33 for example)
    # This is the perihelion velocity.  For a = 1, e = 0, v = 2 pi
    #y[3] = -math.sqrt( (G*M_sun/a_E) * (1.0 + ecc_E) / (1.0 - ecc_E))


    # ================================================================
    # Mars initialization
    # ================================================================

    # our coordinate system puts the primary foci at the origin.

    # set the initial conditions.  The initial position is perihelion
    #y[4] = a_M*(1.0 - ecc_M)  # x 
    #y[5] = 0.0                 # y

    # at perihelion, all the veloicity is in the y-direction.
    #y[6] = 0.0      # v_x

    # v_y^2 = (GM/a) (1+e)/(1-e)  (see C&O Eq. 2.33 for example)
    # This is the perihelion velocity.  For a = 1, e = 0, v = 2 pi
    #y[7] = -math.sqrt( (G*M_sun/a_M) * (1.0 + ecc_M) / (1.0 - ecc_M))


    # These initial conditions were found by putting Mars and Earth
    # both at their perihelion, but with their velocities in the
    # opposite direction (i.e. we want them to go backwards).  This
    # configuration is opposition, and is setup by using the commented
    # out initial conditions above.  We then integrated backwards for 1/4
    # year (91 steps) to get these starting coordinates (note, we reverse
    # the direction of the velocity to get it orbitting in the forward
    # direction.)
    
    # Earth
    y[0] = -0.04631900088483218
    y[1] = -0.9994219951994862
    y[2] = 6.277324691390798
    y[3] = -0.185920887199495

    # Mars
    y[4] = 0.7856599524256417
    y[5] = -1.203323492875661
    y[6] = 4.280834571016523
    y[7] = 3.272064392180777
    

    # set the initial timestep
    dt = 1.0/nsteps_year

    # compute the total number of steps needed 
    nsteps = int(nyears*nsteps_year)

    print "nstep = ", nsteps
    

    # allocate storage to save the positions for all steps
    x_pos_E = numpy.zeros(nsteps, numpy.float64)
    y_pos_E = numpy.zeros(nsteps, numpy.float64)

    x_pos_M = numpy.zeros(nsteps, numpy.float64)
    y_pos_M = numpy.zeros(nsteps, numpy.float64)

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
        x_pos_E[n] = y[0]
        y_pos_E[n] = y[1]

        x_pos_M[n] = y[4]
        y_pos_M[n] = y[5]

        time[n]  = t

        print "step = ", n
        print "Earth: %20.16g %20.16g" % (y[0], y[1])
        print "Earth: %20.16g %20.16g" % (y[2], y[3])
        print " "
        print "Mars:  %20.16g %20.16g" % (y[4], y[5])
        print "Mars:  %20.16g %20.16g" % (y[6], y[7])

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
        pylab.plot(x_pos_E,y_pos_E,color="b")
        pylab.scatter([x_pos_E[n]],[y_pos_E[n]],s=100,color="b")

        # plot planet B
        pylab.plot(x_pos_M,y_pos_M,color="r")
        pylab.scatter([x_pos_M[n]],[y_pos_M[n]],s=100,color="r")

        # draw a line connecting Earth and Mars and extending a bit
        # further out
        slope = (y_pos_M[n] - y_pos_E[n])/(x_pos_M[n] - x_pos_E[n])
        xpt = 3.5
        ypt = y_pos_E[n] + slope*(xpt - x_pos_E[n])

        xline = numpy.zeros(2,numpy.float64)
        yline = numpy.zeros(2,numpy.float64)

        xline[0] = x_pos_E[n]
        xline[1] = xpt
        yline[0] = y_pos_E[n]
        yline[1] = ypt

        pylab.plot(xline,yline,"b--")

        # draw some random background stars
        pylab.scatter([3.2],[ 2.1],s=200,marker=(5,1),color="c")        
        pylab.scatter([3.6],[ 1.0],s=200,marker=(5,1),color="c")        
        pylab.scatter([3.4],[-0.4],s=200,marker=(5,1),color="c")        
        pylab.scatter([3.8],[-0.9],s=200,marker=(5,1),color="c")        
        pylab.scatter([3.1],[-1.3],s=200,marker=(5,1),color="c")        

        pylab.axis([-2.,4.,-2.5,2.5])
        pylab.axis("off")

        f = pylab.gcf()
        f.set_size_inches(6.0,5.0)

        pylab.xlabel("AU")
        pylab.ylabel("AU")
        pylab.text(-1.5,-2.0, "time = %6.3f yr" % time[n])
        pylab.title("Retrograde Mars")

        outfile = "frame_retrograde_%04d.png" % n
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


    
        
