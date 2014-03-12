import math
import numpy
import pylab

# compute the orbits of a planet and the parent star. 
#
# Here we put the center of mass at the origin.  We exaggerate the
# the mass of the planet to make the star's motion more pronounced.

# M. Zingale (2008-11-30)

# we work in CGS units
G = 6.67428e-8        # cm^3 g^{-1} s^{-2}
M_sun = 1.98892e33    # g
AU = 1.49598e13       # cm
year = 3.1557e7

# a simple class to serve as a container for the orbital information
# for the planet and star
class solarsystem:
    
    def __init__ (self, M_star=1, M_p=0):

        self.npts = -1
        self.maxpoints = 2000

        # star's properties
        self.M_star = M_star
        self.x_star = numpy.zeros(self.maxpoints)
        self.y_star = numpy.zeros(self.maxpoints)
        self.vx_star = numpy.zeros(self.maxpoints)
        self.vy_star = numpy.zeros(self.maxpoints)

        # planet's properties
        self.M_p = M_p
        self.x_p = numpy.zeros(self.maxpoints)
        self.y_p = numpy.zeros(self.maxpoints)
        self.vx_p = numpy.zeros(self.maxpoints)
        self.vy_p = numpy.zeros(self.maxpoints)

        self.t = numpy.zeros(self.maxpoints)



def radial_velocity():

    # set the masses
    M_star = M_sun           # star's mass
    M_p    = 0.15*M_sun       # planet's mass

    # set the semi-major axis of the planet (and derive that of the star)
    a_p = 1.0*AU
    a_star = (M_p/M_star)*a_p  # M_p a_p = -M_star a_star (center of mass)


    # create the solar system container
    ss = solarsystem(M_star = M_star, M_p = M_p)


    # set the initial position of the planet
    
    # we are going to put the center of mass at the origin and the planet
    # initially on the +x axis and the star initially on the -x axis
    x_p_init = a_p
    y_p_init = 0.0

    x_star_init = -a_star
    y_star_init = 0.0


    # Kepler's laws should tell us the orbital period
    # P^2 = 4 pi^2 (a_p + a_star)^3 / (G (M_p + M_star))
    # (the semi-major axis of the planet is just |x_p_init|, assuming 
    #  circular orbits; the semi-major axis of the star is |x_star_init|)
    period = math.sqrt(4*math.pi**2*(a_p + a_star)**3/(G*(M_p + M_star)))

    print "period = ", period/year


    # compute the circular velocities
    vx_p_init = 0.0
    vy_p_init = 2*math.pi*a_p/period

    vx_star_init = 0.0
    vy_star_init = -2*math.pi*a_star/period   

    print "star init: ", vx_star_init
    
    # set the timestep in terms of the orbital period
    dt = period/360.0        
    tmax = 2.0*period  # maximum integration time

    integrate_system(ss, 
                     x_p_init, y_p_init, vx_p_init, vy_p_init,
                     x_star_init, y_star_init, vx_star_init, vy_star_init,
                     dt,tmax)


    # ================================================================
    # plotting
    # ================================================================

    # turn on interactive mode 
    #pylab.ion()


    iframe = 0

    # first plot a legend intro sequence
    while (iframe < 100):

        pylab.clf()

        pylab.subplot(121)

        pylab.subplots_adjust(left=0.1,right=0.95,bottom=0.15,top=0.85,
                              wspace=0.45,hspace=0.0)

        a = pylab.gca()
        a.set_aspect("equal", "datalim")
        pylab.axis("off")

        pylab.scatter([0],[0],s=150,marker="x",color="k")


        # plot the planet's orbit and position
        pylab.plot(ss.x_p[0:ss.npts],ss.y_p[0:ss.npts], color="g")
        pylab.scatter([ss.x_p[0]],[ss.y_p[0]],s=100,color="g")

        # plot the star's orbit and position
        pylab.plot(ss.x_star[0:ss.npts],ss.y_star[0:ss.npts], color="r")
        pylab.scatter([ss.x_star[0]],[ss.y_star[0]],s=200,color="r")

        # plot a reference line
        pylab.plot([0,1*AU],[-1.2*a_p,-1.2*a_p], color="k")
        pylab.text(0.5*AU, -1.4*a_p, "1 AU", horizontalalignment='center')

        pylab.arrow(0.0, 1.6*a_p, 0.0, 0.2*a_p, width=0.05*a_p, 
                    facecolor="b",edgecolor="k",alpha=1.0,clip_on=0)
        pylab.text(0.0, 1.4*a_p, "observer", horizontalalignment='center')

        pylab.text(-1.8*a_p, -1.6*a_p, "time = %6.3f yr" % (ss.t[0]/year))

        pylab.axis([-1.5*a_p,1.5*a_p,-1.5*a_p,1.5*a_p])


        pylab.subplot(122)

        
        pylab.scatter([0.1],[0.8], s=200, color="r")
        pylab.text(0.2, 0.8, "star", color="r")

        pylab.scatter([0.1],[0.6], s=100, color="g")
        pylab.text(0.2, 0.6, "planet", color="g")

        pylab.text(0.2, 0.5, "(note: planet's mass exaggerated", size="small",color="k")
        pylab.text(0.2, 0.45, " to enhance effect)", size="small",color="k")


        pylab.scatter([0.1], [0.3], s=150,marker="x",color="k")
        pylab.text(0.2, 0.3, "center of mass", color="k")

        pylab.axis([0,1,0,1])
        pylab.axis("off")

        f = pylab.gcf()
        f.set_size_inches(8.0,4.0)

        outfile = "radial_velocity_%04d.png" % iframe
        pylab.savefig(outfile)

        iframe += 1


    # now plot the actual data
    n = 0
    while (n < ss.npts):

        pylab.clf()

        pylab.subplot(121)

        pylab.subplots_adjust(left=0.1,right=0.95,bottom=0.15,top=0.85,
                              wspace=0.45,hspace=0.0)

        a = pylab.gca()
        a.set_aspect("equal", "datalim")
        pylab.axis("off")

        pylab.scatter([0],[0],s=150,marker="x",color="k")


        # plot the planet's orbit and position
        pylab.plot(ss.x_p[0:ss.npts],ss.y_p[0:ss.npts], color="g")
        pylab.scatter([ss.x_p[n]],[ss.y_p[n]],s=100,color="g")

        # plot the star's orbit and position
        pylab.plot(ss.x_star[0:ss.npts],ss.y_star[0:ss.npts], color="r")
        pylab.scatter([ss.x_star[n]],[ss.y_star[n]],s=200,color="r")

        # plot a reference line
        pylab.plot([0,1*AU],[-1.2*a_p,-1.2*a_p], color="k")
        pylab.text(0.5*AU, -1.4*a_p, "1 AU", horizontalalignment='center')

        pylab.arrow(0.0, 1.6*a_p, 0.0, 0.2*a_p, width=0.05*a_p, 
                    facecolor="b",edgecolor="k", alpha=1.0, clip_on=0)
        pylab.text(0.0, 1.4*a_p, "observer", horizontalalignment='center')

        pylab.text(-1.8*a_p, -1.6*a_p, "time = %6.3f yr" % (ss.t[n]/year))

        pylab.axis([-1.5*a_p,1.5*a_p,-1.5*a_p,1.5*a_p])


        pylab.subplot(122)

        # plot km/s vs years
        # note: radial velocity convention is that if it is moving 
        # toward us, then the velocity is negative.  Since the observer
        # is at +Y, we want to plot -vy as the radial velocity.
        pylab.plot(ss.t[0:ss.npts]/period,-ss.vy_star[0:ss.npts]/1.e5, color="r")
        pylab.scatter([ss.t[n]/period],[-ss.vy_star[n]/1.e5], color="r", s=50)

        # plot a reference line at vy = 0
        pylab.plot([0.0,tmax/period], [0.0,0.0], "k--")


        pylab.axis([0.0, tmax/period,
                    int(numpy.min(ss.vy_star))/1.e5 - 1,
                    int(numpy.max(ss.vy_star))/1.e5 + 1])


        pylab.xlabel("t/period")
        pylab.ylabel("radial velocity [km/s]")

        f = pylab.gcf()
        f.set_size_inches(8.0,4.0)
        
        outfile = "radial_velocity_%04d.png" % iframe
        pylab.savefig(outfile)

        iframe += 1
        n += 1



def integrate_system(ss,
                     x_p_init, y_p_init, vx_p_init, vy_p_init, 
                     x_star_init, y_star_init, vx_star_init, vy_star_init,
                     dt, tmax):


    # allocate storage for R-K intermediate results
    # y[0:3] will hold the star info, y[4:7] will hold the planet info
    k1 = numpy.zeros(8, numpy.float64)
    k2 = numpy.zeros(8, numpy.float64)
    k3 = numpy.zeros(8, numpy.float64)
    k4 = numpy.zeros(8, numpy.float64)

    y = numpy.zeros(8, numpy.float64)
    f = numpy.zeros(8, numpy.float64)



    t = 0.0

    # initial conditions

    # star
    y[0] = x_star_init  # initial star x position
    y[1] = y_star_init  # initial star y position

    y[2] = vx_star_init # initial star x-velocity
    y[3] = vy_star_init # initial star y-velocity

    # planet
    y[4] = x_p_init  # initial planet x position
    y[5] = y_p_init  # initial planet y position

    y[6] = vx_p_init # initial planet x-velocity
    y[7] = vy_p_init # initial planet y-velocity

    ss.x_star[0] = y[0]
    ss.y_star[0] = y[1]

    ss.vx_star[0] = y[2]
    ss.vy_star[0] = y[3]

    ss.x_p[0] = y[4]
    ss.y_p[0] = y[5]

    ss.vx_p[0] = y[6]
    ss.vy_p[0] = y[7]

    ss.t[0] = t

    n = 1
    while (n < ss.maxpoints and t < tmax):
        f = rhs(t, y, ss.M_star, ss.M_p)
        k1[:] = dt*f[:]

        f = rhs(t+0.5*dt, y[:]+0.5*k1[:], ss.M_star, ss.M_p)
        k2[:] = dt*f[:]

        f = rhs(t+0.5*dt, y[:]+0.5*k2[:], ss.M_star, ss.M_p)
        k3[:] = dt*f[:]

        f = rhs(t+dt, y[:]+k3[:], ss.M_star, ss.M_p)
        k4[:] = dt*f[:]

        y[:] += (1.0/6.0)*(k1[:] + 2.0*k2[:] + 2.0*k3[:] + k4[:])

        t = t + dt

        ss.x_star[n] = y[0]
        ss.y_star[n] = y[1]

        ss.vx_star[n] = y[2]
        ss.vy_star[n] = y[3]

        ss.x_p[n] = y[4]
        ss.y_p[n] = y[5]

        ss.vx_p[n] = y[6]
        ss.vy_p[n] = y[7]

        ss.t[n] = t

        n += 1

        #print t, ss.x_p[n], ss.y_p[n]
    
    ss.npts = n



def rhs(t,y,M_star,M_p):

    f = numpy.zeros(8, numpy.float64)

    # y[0] = x_star, y[1] = y_star, y[2] = vx_star, y[3] = vy_star
    # y[4] = x_p,    y[5] = y_p,    y[6] = vx_p,    y[7] = vy_p

    # unpack
    x_star = y[0]
    y_star = y[1]

    vx_star = y[2]
    vy_star = y[3]

    x_p = y[4]
    y_p = y[5]

    vx_p = y[6]
    vy_p = y[7]


    # distance between star and planet
    r = numpy.sqrt((x_p - x_star)**2 + (y_p - y_star)**2)


    f[0] = vx_star  # d(x_star) / dt
    f[1] = vy_star  # d(y_star) / dt

    f[2] = -G*M_p*(x_star - x_p)/r**3  # d(vx_star) / dt
    f[3] = -G*M_p*(y_star - y_p)/r**3  # d(vy_star) / dt

    f[4] = vx_p  # d(x_p) / dt
    f[5] = vy_p  # d(y_p) / dt

    f[6] = -G*M_star*(x_p - x_star)/r**3  # d(vx_p) / dt
    f[7] = -G*M_star*(y_p - y_star)/r**3  # d(vy_p) / dt


    return f
    

    
if __name__== "__main__":
    radial_velocity()


    
        
