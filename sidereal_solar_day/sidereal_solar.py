import math
import numpy
import pylab

# consider a planet whose sidereal rotation period = 1/2 its orbital
# period.  Its solar day would then be 1 orbital period.  This
# animation demonstrates that.


# M. Zingale (2010-10-07)

# we work in CGS units
G = 6.67384e-11      # m^3 kg^{-1} s^{-2}
M_sun = 1.98892e30      # kg
AU = 1.49598261e11         # m
year = 3.155760e7      # s
day = 24*3600           # s


#-----------------------------------------------------------------------------
# a simple class to serve to compute a trajectory of a planet around
# the Sun -- this holds both the orbit information and contains a 
# method to integrate the orbit.
class trajectory:
    
    def __init__ (self):

        self.npts = -1
        self.maxpoints = 2000
        self.x  = numpy.zeros(self.maxpoints)
        self.y  = numpy.zeros(self.maxpoints)
        self.vx = numpy.zeros(self.maxpoints)
        self.vy = numpy.zeros(self.maxpoints)
        self.t  = numpy.zeros(self.maxpoints)

    def integrate(self,x_init,y_init,vx_init,vy_init,dt,tmax):

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
        y[:] = numpy.array([x_init, y_init, vx_init, vy_init])

        # store the initial conditions
        self.x[0] = y[0]
        self.y[0] = y[1]

        self.vx[0] = y[2]
        self.vy[0] = y[3]

        self.t[0] = t

        n = 1
        while (n < self.maxpoints and t < tmax):

            k1[:] = dt*self.rhs(t, y)
            k2[:] = dt*self.rhs(t+0.5*dt, y[:]+0.5*k1[:])
            k3[:] = dt*self.rhs(t+0.5*dt, y[:]+0.5*k2[:])
            k4[:] = dt*self.rhs(t+dt, y[:]+k3[:])

            y[:] += (1.0/6.0)*(k1[:] + 2.0*k2[:] + 2.0*k3[:] + k4[:])

            t = t + dt

            self.x[n]  = y[0]
            self.y[n]  = y[1]
            self.vx[n] = y[2]
            self.vy[n] = y[3]
            self.t[n]  = t

            n += 1

    
        self.npts = n

    def rhs(self,t,y):

        f = numpy.zeros(4, numpy.float64)

        # y[0] = x, y[1] = y, y[2] = v_x, y[3] = v_y
        f[0] = y[2]
        f[1] = y[3]
        r = numpy.sqrt(y[0]**2 + y[1]**2)
        f[2] = -G*M_sun*y[0]/r**3
        f[3] = -G*M_sun*y[1]/r**3

        return f


#-----------------------------------------------------------------------------
def sidereal():

    # set the semi-major axis and eccentricity
    a = 1.0*AU
    e = 0.0

    # compute the period of the orbit from Kepler's law
    P_orbital = math.sqrt(4*math.pi*math.pi*a**3/(G*M_sun))


    # set the rotation period -- this is the sidereal day
    #P_rotation = 0.99726968*day    # Earth
    P_rotation = 0.05*P_orbital

    omega = 2*math.pi/P_rotation

    # set the initial coordinates -- perihelion
    x_init = 0.0
    y_init = -a*(1.0 - e)

    # set the initial velocity
    vx_init = math.sqrt( (G*M_sun/a) * (1 + e) / (1 - e))
    vy_init = 0.0


    # We will maintain two trajectories: full_orbit is the full orbit
    # of the planet.  orbit is just a small segment carrying 2 
    # sidereal days
    orbit = trajectory()
    full_orbit = trajectory()

    # compute the length of the solar day
    P_solar = P_rotation/(1.0 - P_rotation/P_orbital)

    print "period = ", P_orbital/year

    print "sidereal day = ", P_rotation/day
    print "solar day    = ", P_solar/day


    # evolve for 2 sidereral days with lots of frames to get a smooth
    # animation
    tmax = 2.0*P_rotation
    dt = tmax/1000.0
    orbit.integrate(x_init,y_init,vx_init,vy_init,dt,tmax)


    # evolve the full orbit for plotting purposes
    tmax = P_orbital
    dt_f = tmax/720.0
    full_orbit.integrate(x_init,y_init,vx_init,vy_init,dt_f,tmax)    


    # ================================================================
    # plotting
    # ================================================================

    # plot the orbit
    iframe = 0

    # v1
    n = 0
    iframe = 0

    while (n < orbit.npts):

        pylab.clf()


        # plot the foci
        pylab.scatter([0],[0],s=1600,marker=(20,1),color="k")
        pylab.scatter([0],[0],s=1500,marker=(20,1),color="#FFFF00")

        # plot the orbit
        pylab.plot(full_orbit.x[0:full_orbit.npts],
                   full_orbit.y[0:full_orbit.npts],color="r")


        # plot planet 
        theta = numpy.arange(180)
        r = 0.05*AU  # exaggerate the planet's size
        x_surface = orbit.x[n] + r*numpy.cos(theta)
        y_surface = orbit.y[n] + r*numpy.sin(theta)
        pylab.fill(x_surface,y_surface,"b", edgecolor="b",zorder=1000)

        # plot a point on the planet's surface
        xpt = orbit.x[n] + 1.2*r*numpy.cos(omega*orbit.t[n]+math.pi/2.0)
        ypt = orbit.y[n] + 1.2*r*numpy.sin(omega*orbit.t[n]+math.pi/2.0)
        pylab.scatter([xpt],[ypt],s=15,color="k")


        pylab.axis([-1.2*AU,1.2*AU,-1.2*AU,1.2*AU])

        f = pylab.gcf()
        f.set_size_inches(7.2, 7.2)

        pylab.title("Sidereal vs. Solar Day")

        pylab.axis("off")

        ax = pylab.gca()
        ax.set_aspect("equal", "datalim")


        pylab.text(-1.0*AU, -1.2*AU, "Note: length of day exaggerated", 
                    fontsize=9, color="0.5")


        if (orbit.t[n] >= P_rotation - 0.5*dt and
            orbit.t[n] <= P_rotation + 0.5*dt):
            n_sidereal = n


        # special lines
        if (orbit.t[n] >= P_rotation - 0.5*dt and 
            orbit.t[n] <= P_solar + 0.5*dt):
            pylab.plot([orbit.x[n_sidereal], orbit.x[n_sidereal]],
                       [orbit.y[n_sidereal], 0], linestyle=":", color="0.5")


        # special notes
        if (orbit.t[n] < 1.5*dt):
            pylab.text(-1.0*AU, -1.15*AU, "Noon (Sun is on the meridian)")

            d = 0
            while (d < 100):
                outfile = "rotation_%04d.png" % iframe
                pylab.savefig(outfile)            
                iframe += 1
                d += 1

        elif (orbit.t[n] >= P_rotation - 0.5*dt and 
              orbit.t[n] <= P_rotation + 0.5*dt):
            pylab.text(-1.0*AU, -1.15*AU, "1 Sidereal period (Earth rotated 360 degrees)")

            d = 0
            while (d < 100):
                outfile = "rotation_%04d.png" % iframe
                pylab.savefig(outfile)            
                iframe += 1
                d += 1

        elif (orbit.t[n] >= P_solar - 0.5*dt and 
              orbit.t[n] <= P_solar + 0.5*dt):
            pylab.text(-1.0*AU, -1.15*AU, "1 Solar period (Sun is back on the meridian)")

            d = 0
            while (d < 100):
                outfile = "rotation_%04d.png" % iframe
                pylab.savefig(outfile)            
                iframe += 1
                d += 1


        outfile = "rotation_%04d.png" % iframe
        pylab.savefig(outfile)

        iframe += 1

        n += 1
    

    
if __name__== "__main__":
    sidereal()


    
        
