import math
import numpy
import pylab

# Integrate the orbit of Mercury around the Sun and show its rotation 
# by drawing a point on the surface.


# M. Zingale (2011-08-31)

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




def orbitalenergy():

    # set the semi-major axis and eccentricity
    a = 0.387*AU
    e = 0.2056

    # compute the period of the orbit from Kepler's law and make
    # the timestep by 1/720th of a period
    P_orbital = math.sqrt(4*math.pi*math.pi*a**3/(G*M_sun))


    # set the rotation period
    P_rotation = (2./3.)*P_orbital

    omega = 2*math.pi/P_rotation

    # set the initial coordinates -- perihelion
    x_init = a*(1.0 - e)
    y_init = 0.0

    # set the initial velocity
    vx_init = 0.0
    vy_init = math.sqrt( (G*M_sun/a) * (1 + e) / (1 - e))


    # the circular orbit will be our baseline
    orbit = trajectory()


    print "period = ", P_orbital/year

    dt = P_orbital/360.0
    tmax = 4*P_orbital

    integrate_projectile(orbit,x_init,y_init,vx_init,vy_init,dt,tmax)




    # ================================================================
    # plotting
    # ================================================================

    # turn on interactive mode 
    #pylab.ion()


    # plot the orbit
    iframe = 0

    # v1
    n = 0
    while (n < orbit.npts):

        pylab.clf()


        # plot the foci
        pylab.scatter([0],[0],s=250,marker=(5,1),color="k")
        pylab.scatter([0],[0],s=200,marker=(5,1),color="y")

        # plot the orbit
        pylab.plot(orbit.x[0:orbit.npts],orbit.y[0:orbit.npts],color="r")
        #pylab.scatter([orbit.x[n]],[orbit.y[n]],s=100,color="r")

        # plot planet 
        theta = numpy.arange(180)
        r = 0.03*AU  # exaggerate the planet's size
        x_surface = orbit.x[n] + r*numpy.cos(theta)
        y_surface = orbit.y[n] + r*numpy.sin(theta)
        pylab.fill(x_surface,y_surface,color="r", alpha=1.0, zorder=1000)

        # plot a point on the planet's surface
        xpt = orbit.x[n] + r*numpy.cos(omega*orbit.t[n]+math.pi)
        ypt = orbit.y[n] + r*numpy.sin(omega*orbit.t[n]+math.pi)
        pylab.scatter([xpt],[ypt],s=25,color="k")


        pylab.axis([-0.75*AU,0.75*AU,-0.75*AU,0.75*AU])

        f = pylab.gcf()
        f.set_size_inches(6.0,6.0)


        pylab.xlabel("x [m]")
        pylab.ylabel("y [m]")


        outfile = "rotation_%04d.png" % n
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


    
        
