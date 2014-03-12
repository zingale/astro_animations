import math
import numpy
import pylab

# compute the orbit of a planet around the Sun using fourth-order accurate
# Runge-Kutta

# we work in units of AU, yr, and M_sun
#  in these units, G = 4 pi^2

pi = math.pi
G = 4.0*pi*pi  # AU^3 M_sun^{-1} yr^{-2}
M_sun = 1.0
ecc = 0.4     # eccentricity
a = 1.0       # semi-major axis


def integrate():

    # allocate storage
    k1 = numpy.zeros(4, numpy.float64)
    k2 = numpy.zeros(4, numpy.float64)
    k3 = numpy.zeros(4, numpy.float64)
    k4 = numpy.zeros(4, numpy.float64)

    # y[0], y[1] are (x,y); y[2], y[3] are (v_x, v_y)
    y = numpy.zeros(4, numpy.float64)
    f = numpy.zeros(4, numpy.float64)

    t = 0.0
    dt = 0.0


    # set the initial conditions.  The initial position is  perihelion
    y[0] = a*(1.0 - ecc)   # x -- this is the semi-major axis
    y[1] = 0.0             # y

    # at perihelion, all the veloicity is in the y-direction.
    y[2] = 0.0      # v_x

    # v_y^2 = (GM/a) (1+e)/(1-e)  (see C&O Eq. 2.33 for example)
    # This is the perihelion velocity.  For a = 1, e = 0, v = 2 pi
    y[3] = math.sqrt( (G*M_sun/a) * (1.0 + ecc) / (1.0 - ecc))

    # set the initial timestep
    nsteps = 365

    dt = 1.0/nsteps

    # compute the number of steps needed for 2 years of integration
    nsteps = 2*nsteps


    # allocate storage to save the positions for all steps
    x_pos = numpy.zeros(nsteps, numpy.float64)
    y_pos = numpy.zeros(nsteps, numpy.float64)
    time  = numpy.zeros(nsteps, numpy.float64)

    #print "# t, x, y, vx, vy"

    #print t, y[0], y[1], y[2], y[3]

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

        x_pos[n] = y[0]
        y_pos[n] = y[1]
        time[n]  = t

        t = t + dt
        n += 1



        #print t, y[0], y[1], y[2], y[3]


    
    # turn on interactive mode 
    pylab.ion()

    # now loop over and plot point by point
    #n = 0
    #while (n < nsteps):
    #    fname = '_tmp%03d.png'%n
    #    pylab.clf()
    #    pylab.axis([-1.5,1.0,-1.25,1.25])
    #    pylab.plot(x_pos,y_pos)
    #    pylab.plot([0,x_pos[n]],[0,y_pos[n]])
    #    pylab.scatter([0],[0])
    #    pylab.scatter([x_pos[n]],[y_pos[n]])
    #    pylab.draw()
    #    n += 1

    # finally loop over and show equal area in equal time
    intervals = 12
    t_interval = (0.5*t/intervals)*numpy.arange(intervals+1)
    mask = 1-numpy.mod(numpy.arange(intervals),3)

    print t_interval
    print mask

    pylab.axis([-1.5,1.0,-1.25,1.25])
    pylab.plot(x_pos,y_pos)

    n = 0
    while (n < nsteps/2):
        fname = '_tmp%03d.png'%n


        # what interval are we in
        i = 0
        intv = -1
        while (i < intervals):
            if (time[n] > t_interval[i] and time[n] <= t_interval[i+1]):
                intv = i
                break
            
            i += 1

            
        if (intv >= 0 and mask[intv] == 1):
            pylab.plot([0,x_pos[n]],[0,y_pos[n]])

        pylab.scatter([0],[0])
        pylab.draw()
        n += 1




def rhs(t,y):

    f = numpy.zeros(4, numpy.float64)

    # y[0] = x, y[1] = y, y[2] = v_x, y[3] = v_y
    f[0] = y[2]
    f[1] = y[3]
    r = numpy.sqrt(y[0]**2 + y[1]**2)
    f[2] = -G*M_sun*y[0]/r**3
    f[3] = -G*M_sun*y[1]/r**3

    return f


    
if __name__== "__main__":
    integrate()


    
        
