import math
import numpy
import pylab

# illustrate Kepler's second law by showing the areas corresponding to
# equal time intervals.

# M. Zingale (2009-03-01)

# we work in MKS units
G = 6.67428e-11      # m^3 kg^{-1} s^{-2}
M_sun = 1.98892e30   # kg
AU = 1.49598e11      # m
year = 3.1557e7      # s


# a simple class to serve as a container for the orbital information
class trajectory:
    
    def __init__ (self):

        self.npts = -1
        self.maxpoints = 2000
        self.x  = numpy.zeros(self.maxpoints)
        self.y  = numpy.zeros(self.maxpoints)
        self.vx = numpy.zeros(self.maxpoints)
        self.vy = numpy.zeros(self.maxpoints)
        self.t  = numpy.zeros(self.maxpoints)


class interval:

    def __init__ (self, start, end):

        self.start = start
        self.end = end


def orbitalenergy():

    # set the semi-major axis and eccentricity
    a = 1.5874*AU
    e = 0.5

    # set the initial coordinates -- perihelion
    x_init = a*(1.0 - e)
    y_init = 0.0

    # set the initial velocity
    vx_init = 0.0
    vy_init = math.sqrt( (G*M_sun/a) * (1 + e) / (1 - e))


    # the circular orbit will be our baseline
    orbit = trajectory()

    # compute the period of the orbit from Kepler's law and make
    # the timestep by 1/720th of a period
    P = math.sqrt(4*math.pi*math.pi*a**3/(G*M_sun))

    print "period = ", P/year

    dt = P/720.0
    tmax = P

    integrate_projectile(orbit,x_init,y_init,vx_init,vy_init,dt,tmax)

    
    # set up some intervals (start time, end time) to shade.  
    # Note, the length of the intervals should be the same
    intervals = []
    intervals.append(interval(0, P/12))
    intervals.append(interval(P/2, 7.0*P/12.0))
    intervals.append(interval(9.5*P/12, 10.5*P/12.0))

    i = 0
    while (i < len(intervals)):

        if (intervals[i].start < 0 or intervals[i].end > P):
            print "ERROR: interval %d is not contained in a single orbit" % \
                (i)

        print "interval %d, dt = %f" % \
            (i, intervals[i].end - intervals[i].start)
        i += 1


    # ================================================================
    # plotting
    # ================================================================

    # plot the orbit
    iframe = 0

    # v1
    n = 0
    while (n < orbit.npts):

        pylab.clf()

        pylab.subplots_adjust(left=0.1,right=0.9,bottom=0.1,top=0.9)

        a = pylab.gca()
        a.set_aspect("equal", "datalim")
        pylab.axis("off")

        # plot the foci
        pylab.scatter([0],[0],s=250,marker=(5,1),color="k")
        pylab.scatter([0],[0],s=200,marker=(5,1),color="y")

        # plot planet 
        pylab.plot(orbit.x[0:orbit.npts],orbit.y[0:orbit.npts],color="r")
        pylab.scatter([orbit.x[n]],[orbit.y[n]],s=100,color="r")
        
        # shaded the intervals
        i = 0
        while (i < len(intervals)):

            x = [0]
            y = [0]

            m = 0
            while (m < orbit.npts):

                if (orbit.t[m] >= intervals[i].start and
                    orbit.t[m] <= min(intervals[i].end,orbit.t[n])):
                    x.append(orbit.x[m])
                    y.append(orbit.y[m])

                elif (orbit.t[m] > intervals[i].end):
                    break

                m += 1

            pylab.fill(x,y,alpha=0.20, facecolor="b")

            i += 1


        pylab.axis([-4*AU,2*AU,-2*AU,2*AU])        

        f = pylab.gcf()
        f.set_size_inches(6.0,4.0)

        pylab.xlabel("x [m]")
        pylab.ylabel("y [m]")

        outfile = "second_law_%04d.png" % n
        pylab.savefig(outfile)
        n += 1




def integrate_projectile(orbit,x_init,y_init,vx_init,vy_init,dt,tmax):

    SMALL = 1.e-16

    # allocate storage for R-K intermediate results
    k1 = numpy.zeros(4, numpy.float64)
    k2 = numpy.zeros(4, numpy.float64)
    k3 = numpy.zeros(4, numpy.float64)
    k4 = numpy.zeros(4, numpy.float64)


    y = numpy.zeros(4, numpy.float64)
    f = numpy.zeros(4, numpy.float64)



    t = 0.0

    # initial conditions
    y[0] = x_init   
    y[1] = y_init   

    y[2] = vx_init   
    y[3] = vy_init     

    # store the initial conditions
    orbit.x[0] = y[0]
    orbit.y[0] = y[1]

    orbit.vx[0] = y[2]
    orbit.vy[0] = y[3]

    orbit.t[0] = t


    n = 1
    while (n < orbit.maxpoints and t < tmax):

        f = rhs(t, y)
        k1[:] = dt*f[:]

        f = rhs(t+0.5*dt, y[:]+0.5*k1[:])
        k2[:] = dt*f[:]

        f = rhs(t+0.5*dt, y[:]+0.5*k2[:])
        k3[:] = dt*f[:]

        f = rhs(t+dt, y[:]+k3[:])
        k4[:] = dt*f[:]

        y[:] += (1.0/6.0)*(k1[:] + 2.0*k2[:] + 2.0*k3[:] + k4[:])

        t = t + dt

        orbit.x[n]  = y[0]
        orbit.y[n]  = y[1]
        orbit.vx[n] = y[2]
        orbit.vy[n] = y[3]
        orbit.t[n]  = t

        n += 1

    
    orbit.npts = n



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
    orbitalenergy()


    
        
