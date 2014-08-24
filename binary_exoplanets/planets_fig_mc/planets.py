import math
import numpy
import pylab

# compute the orbit of a 2 planets around the Sun using fourth-order accurate
# Runge-Kutta
#
# show two planets with the same perihelion, but different eccentricities
#
# M. Zingale (2010-09-16)

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

# we want the perihelion distance of planet B to be the same
# as planet a
r_p_A = a_A*(1.0 - ecc_A)
a_B = r_p_A/(1.0 - ecc_B)


# xmin and xmax
xmin = min(-a_A*(1.0 + ecc_A), -a_B*(1.0 + ecc_B) )
xmax = max( a_A*(1.0 - ecc_A),  a_B*(1.0 - ecc_B) )



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

    # set the initial conditions.  The initial position is perihelion
    y[4] = a_B*(1.0 - ecc_B)   # x 
    y[5] = 0.0                 # y

    # at aphelion, all the veloicity is in the negative y-direction.
    y[6] = 0.0      # v_x

    # v_y^2 = (GM/a) (1+e)/(1-e)  (see C&O Eq. 2.33 for example)
    # This is the perihelion velocity.  For a = 1, e = 0, v = 2 pi
    y[7] = math.sqrt( (G*M_sun/a_B) * (1.0 + ecc_B) / (1.0 - ecc_B))



    # set the initial timestep
    dt = 1.0/nsteps_year


    # compute the periods
    P_A = math.sqrt( 4.0*math.pi**2 * a_A**3 / (G*M_sun) )
    P_B = math.sqrt( 4.0*math.pi**2 * a_B**3 / (G*M_sun) )


    # compute the total number of steps needed 
    nsteps = int(nyears*nsteps_year)

    # compute the number of steps for one period
    nsteps_A = int(P_A/dt) + 1
    nsteps_B = int(P_B/dt) + 1


    print "period A: ", P_A, nsteps_A
    print "period_B: ", P_B, nsteps_B


    # allocate storage to save the positions for all steps
    x_pos_A = numpy.zeros(nsteps, numpy.float64)
    y_pos_A = numpy.zeros(nsteps, numpy.float64)

    x_pos_B = numpy.zeros(nsteps, numpy.float64)
    y_pos_B = numpy.zeros(nsteps, numpy.float64)

    time  = numpy.zeros(nsteps, numpy.float64)


    # store the initial conditions
    x_pos_A[0] = y[0]
    y_pos_A[0] = y[1]
    
    x_pos_B[0] = y[4]
    y_pos_B[0] = y[5]

    n = 1



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

    pylab.clf()

    # plot the foci
    pylab.scatter([0],[0],s=250,marker=(5,1),color="k")
    pylab.scatter([0],[0],s=200,marker=(5,1),color="y")

    # plot planet A orbit
    pylab.plot(x_pos_A,y_pos_A,color="r")

    # plot planet B orbit
    pylab.plot(x_pos_B,y_pos_B,color="b",linestyle="--")

    # plot planet B points
    pylab.scatter([x_pos_B[0]],[y_pos_B[0]],s=100,color="b",marker="h")
    pylab.text(x_pos_B[0]*1.15,y_pos_B[0], "B1", color="b")

    pylab.scatter([x_pos_B[nsteps_B/2]],[y_pos_B[nsteps_B/2]],s=100,color="b",marker="h")
    pylab.text(x_pos_B[nsteps_B/2]*1.15,y_pos_B[nsteps_B/2], "B2", color="b")

    pylab.scatter([x_pos_B[3*nsteps_B/4]],[y_pos_B[3*nsteps_B/4]],s=100,color="b",marker="h")
    pylab.text(x_pos_B[3*nsteps_B/4],y_pos_B[3*nsteps_B/4]*1.15, "B3", color="b")

    # find the geometric 1/2 point between aphelion and perihelion
    x_peri = x_pos_B[0]
    y_peri = y_pos_B[0]

    x_aph = x_pos_B[nsteps_B/2]
    y_aph = y_pos_B[nsteps_B/2]

    n = nsteps_B/2
    

    while (n < nsteps_B):
        
        d_from_peri = math.sqrt( (x_pos_B[n] - x_peri)**2 + (y_pos_B[n] - y_peri)**2 )
        d_from_aph  = math.sqrt( (x_pos_B[n] - x_aph)**2  + (y_pos_B[n] - y_aph)**2 )
        
        # we start out closer to aphelion.  As soon as that changes,
        # we are at the 1/2-way point
        if (d_from_peri <= d_from_aph):
            n_half_B = n
            break
        
        n += 1

    pylab.scatter([x_pos_B[n_half_B]],[y_pos_B[n_half_B]],s=100,color="b",marker="h")
    pylab.text(x_pos_B[n_half_B],y_pos_B[n_half_B]*1.15, "B4", color="b")

    # last point equidistance spacing after the geometric 1/2 point
    n_next = (n_half_B - 3*nsteps_B/4) + n_half_B

    pylab.scatter([x_pos_B[n_next]],[y_pos_B[n_next]],s=100,color="b",marker="h")
    pylab.text(x_pos_B[n_next],y_pos_B[n_next]*1.15, "B5", color="b")    


    # plot planet A points
    pylab.scatter([x_pos_A[0]],[y_pos_A[0]],s=100,color="r",marker="^")
    pylab.text(x_pos_A[0]*0.7,y_pos_A[0], "A1", color="r")



    # draw the center of the ellipse
    pylab.scatter([-a_B*ecc_B], [0], s=100, color='b', marker='x')

    # draw the semi-major and semi-minor axes
    pylab.plot([-a_B*ecc_B,-a_B*ecc_B], [-a_B*math.sqrt(1.0 - ecc_B**2),a_B*math.sqrt(1.0 - ecc_B**2)], color='b', linestyle=":")
    pylab.plot([-a_B*(1.0 + ecc_B), a_B*(1.0 - ecc_B)], [0,0], color='b', linestyle=":")

    ax = pylab.gca()
    ax.set_aspect("equal", "datalim")
    pylab.axis("off")

    pylab.axis([1.1*xmin,1.5*xmax,-0.5*(xmax - xmin), 0.5*(xmax - xmin)])
    
    f = pylab.gcf()
    f.set_size_inches(6.0,6.0)

    pylab.xlabel("AU")
    pylab.ylabel("AU")
    
    outfile = "planets.eps"
    pylab.savefig(outfile)



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


    
        
