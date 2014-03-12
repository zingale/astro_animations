import math
import numpy
import pylab

# compute the orbits of stars in a binary system.
#
# Here we put the star 1 at the origin -- this shows what the orbit of
# the second star looks like from the point of view of the first.
#
# This version allows for elliptical orbits with some arbitrary
# orientation wrt to the observer (although, still face-on)
#
# The plotting assumes that M_2 < M_1
#  
# M. Zingale (2009-02-12)

# we work in CGS units
G = 6.67428e-8        # cm^3 g^{-1} s^{-2}
M_sun = 1.98892e33    # g
AU = 1.49598e13       # cm
year = 3.1556926e7

# a simple class to serve as a container for the orbital information
# for the two stars
class solarsystem:
    
    def __init__ (self, M_star1=1, M_star2=0.1):

        self.npts = -1
        self.maxpoints = 2000

        # star1 properties
        self.M_star1 = M_star1
        self.x_star1 = numpy.zeros(self.maxpoints)
        self.y_star1 = numpy.zeros(self.maxpoints)
        self.vx_star1 = numpy.zeros(self.maxpoints)
        self.vy_star1 = numpy.zeros(self.maxpoints)

        # planet's properties
        self.M_star2 = M_star2
        self.x_star2 = numpy.zeros(self.maxpoints)
        self.y_star2 = numpy.zeros(self.maxpoints)
        self.vx_star2 = numpy.zeros(self.maxpoints)
        self.vy_star2 = numpy.zeros(self.maxpoints)

        self.t = numpy.zeros(self.maxpoints)


def find_scinotat(number):

    b = int(math.log10(math.fabs(number)))
    a = number/10**b
    
    return a, b


def radial_velocity():

    # set the masses
    M_star1 = 4.0*M_sun      # star 1's mass
    M_star2 = M_sun      # star 2's mass

    # set the semi-major axis of the star 2 (and derive that of star 1)
    # M_star2 a_star2 = -M_star1 a_star1 (center of mass)
    a_star2 = 10.0*AU
    a_star1 = (M_star2/M_star1)*a_star2  

    # set the eccentricity
    ecc = 0.4

    # set the angle to rotate the semi-major axis wrt the observer
    theta = math.pi/6.0

    # display additional information
    annotate = 1


    # create the solar system container
    ss = solarsystem(M_star1 = M_star1, M_star2 = M_star2)


    # set the initial position of the planet -- perihelion
    
    # we are going to put the center of mass at the origin and star 2
    # initially on the +x axis and the star 1 initially on the -x axis
    x_star1_init = -a_star1*(1.0 - ecc)*math.cos(theta)
    y_star1_init = -a_star1*(1.0 - ecc)*math.sin(theta)

    x_star2_init = a_star2*(1.0 - ecc)*math.cos(theta)
    y_star2_init = a_star2*(1.0 - ecc)*math.sin(theta)


    # Kepler's laws should tell us the orbital period
    # P^2 = 4 pi^2 (a_star1 + a_star2)^3 / (G (M_star1 + M_star2))
    period = math.sqrt(4*math.pi**2*(a_star1 + a_star2)**3/(G*(M_star1 + M_star2)))

    print "period = ", period/year


    # compute the velocities.

    # first compute the velocity of the reduced mass at perihelion
    # (C&O Eq. 2.33)
    v_mu = math.sqrt( (G*(M_star1 + M_star2)/(a_star1 + a_star2)) *
                      (1.0 + ecc)/(1.0 - ecc) )

    # then v_star2 = (mu/m_star2)*v_mu
    vx_star2_init = -(M_star1/(M_star1 + M_star2))*v_mu*math.sin(theta)
    vy_star2_init = (M_star1/(M_star1 + M_star2))*v_mu*math.cos(theta)

    # then v_star1 = (mu/m_star1)*v_mu
    vx_star1_init = (M_star2/(M_star1 + M_star2))*v_mu*math.sin(theta)
    vy_star1_init = -(M_star2/(M_star1 + M_star2))*v_mu*math.cos(theta)

    
    # set the timestep in terms of the orbital period
    dt = period/360.0        
    tmax = 2.0*period  # maximum integration time

    integrate_system(ss, 
                     x_star1_init, y_star1_init, vx_star1_init, vy_star1_init,
                     x_star2_init, y_star2_init, vx_star2_init, vy_star2_init,
                     dt,tmax)


    # ================================================================
    # plotting
    # ================================================================

    # turn on interactive mode 
    #pylab.ion()


    iframe = 0

    n = 0
    while (n < ss.npts):
        
        pylab.clf()

        pylab.rc('text', usetex=True)
        pylab.rc('font', family='serif')

        pylab.subplots_adjust(left=0.1,right=0.9,bottom=0.1,top=0.9)

        a = pylab.gca()
        a.set_aspect("equal", "datalim")
        pylab.axis("off")


        # offsets -- we want star 1 at the origin
        xoffset = ss.x_star1[n]
        yoffset = ss.y_star1[n]

        pylab.scatter([0-xoffset],[0-yoffset],s=150,marker="x",color="k")

        # if e = 0 and M_star1 = M_star2, then the orbits lie on top of one
        # another, so plot only a single orbital line.

        # plot star 1's position
        symsize = 200
        pylab.scatter([ss.x_star1[n]-xoffset],[ss.y_star1[n]-yoffset],s=symsize,color="r")

        # plot star 2's orbit and position
        symsize = 200*(M_star2/M_star1)
        pylab.plot(ss.x_star2[0:ss.npts]-ss.x_star1[0:ss.npts],
                   ss.y_star2[0:ss.npts]-ss.y_star1[0:ss.npts], color="g")

        pylab.scatter([ss.x_star2[n]-xoffset],[ss.y_star2[n]-yoffset],s=symsize,color="g")


        if (annotate):
            # plot a reference line
            pylab.plot([0,1*AU],[-1.6*a_star2,-1.6*a_star2], color="k")
            pylab.text(0.5*AU, -1.8*a_star2, "1 AU", 
                       horizontalalignment='center')

            # display time
            pylab.text(-1.8*a_star2, -1.7*a_star2, 
                        "time = %6.3f yr" % (ss.t[n]/year))

        # display information about stars
        pylab.text(-1.8*a_star2, 0.9*a_star2, 
                    r"mass ratio: %3.2f" % (ss.M_star1/ss.M_star2), color="k")
        pylab.text(-1.8*a_star2, 0.7*a_star2, 
                    r"eccentricity: %3.2f" % (ecc), color="k")



        pylab.axis([-1.8*a_star2,1.0*a_star2,-1.8*a_star2,1.0*a_star2])

        f = pylab.gcf()
        f.set_size_inches(6.0,6.0)
  
        outfile = "binary_star_%04d.png" % iframe
        pylab.savefig(outfile)


        iframe += 1
        n += 1



def integrate_system(ss,
                     x_star1_init, y_star1_init, vx_star1_init, vy_star1_init,
                     x_star2_init, y_star2_init, vx_star2_init, vy_star2_init, 
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
    y[0] = x_star1_init  # initial star x position
    y[1] = y_star1_init  # initial star y position

    y[2] = vx_star1_init # initial star x-velocity
    y[3] = vy_star1_init # initial star y-velocity

    # planet
    y[4] = x_star2_init  # initial planet x position
    y[5] = y_star2_init  # initial planet y position

    y[6] = vx_star2_init # initial planet x-velocity
    y[7] = vy_star2_init # initial planet y-velocity

    ss.x_star1[0] = y[0]
    ss.y_star1[0] = y[1]

    ss.vx_star1[0] = y[2]
    ss.vy_star1[0] = y[3]

    ss.x_star2[0] = y[4]
    ss.y_star2[0] = y[5]

    ss.vx_star2[0] = y[6]
    ss.vy_star2[0] = y[7]

    ss.t[0] = t

    n = 1
    while (n < ss.maxpoints and t < tmax):

        f = rhs(t, y, ss.M_star1, ss.M_star2)
        k1[:] = dt*f[:]

        f = rhs(t+0.5*dt, y[:]+0.5*k1[:], ss.M_star1, ss.M_star2)
        k2[:] = dt*f[:]

        f = rhs(t+0.5*dt, y[:]+0.5*k2[:], ss.M_star1, ss.M_star2)
        k3[:] = dt*f[:]

        f = rhs(t+dt, y[:]+k3[:], ss.M_star1, ss.M_star2)
        k4[:] = dt*f[:]

        y[:] += (1.0/6.0)*(k1[:] + 2.0*k2[:] + 2.0*k3[:] + k4[:])

        t = t + dt

        ss.x_star1[n] = y[0]
        ss.y_star1[n] = y[1]

        ss.vx_star1[n] = y[2]
        ss.vy_star1[n] = y[3]

        ss.x_star2[n] = y[4]
        ss.y_star2[n] = y[5]

        ss.vx_star2[n] = y[6]
        ss.vy_star2[n] = y[7]

        ss.t[n] = t

        n += 1

        #print t, ss.x_star2[n], ss.y_star2[n]
    
    ss.npts = n



def rhs(t,y,M_star1,M_star2):

    f = numpy.zeros(8, numpy.float64)

    # y[0] = x_star1, y[1] = y_star1, y[2] = vx_star1, y[3] = vy_star1
    # y[4] = x_star2,    y[5] = y_star2,    y[6] = vx_star2,    y[7] = vy_star2

    # unpack
    x_star1 = y[0]
    y_star1 = y[1]

    vx_star1 = y[2]
    vy_star1 = y[3]

    x_star2 = y[4]
    y_star2 = y[5]

    vx_star2 = y[6]
    vy_star2 = y[7]


    # distance between star and planet
    r = numpy.sqrt((x_star2 - x_star1)**2 + (y_star2 - y_star1)**2)


    f[0] = vx_star1  # d(x_star1) / dt
    f[1] = vy_star1  # d(y_star1) / dt

    f[2] = -G*M_star2*(x_star1 - x_star2)/r**3  # d(vx_star1) / dt
    f[3] = -G*M_star2*(y_star1 - y_star2)/r**3  # d(vy_star1) / dt

    f[4] = vx_star2  # d(x_star2) / dt
    f[5] = vy_star2  # d(y_star2) / dt

    f[6] = -G*M_star1*(x_star2 - x_star1)/r**3  # d(vx_star2) / dt
    f[7] = -G*M_star1*(y_star2 - y_star1)/r**3  # d(vy_star2) / dt


    return f
    

    
if __name__== "__main__":
    radial_velocity()


    
        
