import math
import numpy
import pylab

# compute the orbits of eclipsing binary stars, seen nearly edge-on.
#
# Here we put the center of mass at the origin.  
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
deg_to_rad = math.pi/180.0

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

        # star2 properties
        self.M_star2 = M_star2
        self.x_star2 = numpy.zeros(self.maxpoints)
        self.y_star2 = numpy.zeros(self.maxpoints)
        self.vx_star2 = numpy.zeros(self.maxpoints)
        self.vy_star2 = numpy.zeros(self.maxpoints)

        self.t = numpy.zeros(self.maxpoints)



def radial_velocity():

    # set the masses
    M_star1 = 5*M_sun           # star 1's mass
    M_star2 = M_sun           # star 2's mass

    # set the semi-major axis of the star 2 (and derive that of star 1)
    # M_star2 a_star2 = -M_star1 a_star1 (center of mass)
    a_star2 = 1.0*AU
    a_star1 = (M_star2/M_star1)*a_star2  

    # set the eccentricity
    ecc = 0.0

    # set the angle to rotate the semi-major axis wrt the observer
    theta = 0.0

    # set the incliination wrt the observer.
    inc = 88.75    # degrees

    # create the solar system container
    ss = solarsystem(M_star1 = M_star1, M_star2 = M_star2)


    # set the radii and temperatures -- dimensionless -- we are going
    # to exaggerate their scale
    R1 = 20
    R2 = 5

    rad_scal = 0.01*AU


    # set the temperatures
    T1 = 5000 # K
    T2 = 10000  # K


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
    tmax = period  # maximum integration time

    integrate_system(ss, 
                     x_star1_init, y_star1_init, vx_star1_init, vy_star1_init,
                     x_star2_init, y_star2_init, vx_star2_init, vy_star2_init,
                     dt,tmax)


    # apply the projection to account for the inclination wrt the
    # observer
    ss.y_star1[0:ss.npts] = ss.y_star1[0:ss.npts]*math.cos(inc*deg_to_rad)
    ss.y_star2[0:ss.npts] = ss.y_star2[0:ss.npts]*math.cos(inc*deg_to_rad)

    # we will keep one star fixed and draw the relative `orbit' of the
    # other
    x_rel = ss.x_star1 - ss.x_star2
    y_rel = ss.y_star1 - ss.y_star2


    # cut angle -- we don't want to plot the orbit where the stationary
    # star is (origin), so specify the angle wrt the +y axis where we
    # don't plot
    cut_angle = 10*deg_to_rad

    frac = cut_angle/(2.0*math.pi)

    # star 2 starts on the -x axis.  3/4 through the orbit, it will
    # be behind the star 1.  Compute the range of steps to skip plotting
    cut_index1 = 0.75*ss.npts - frac*ss.npts
    cut_index2 = 0.75*ss.npts + frac*ss.npts

    print cut_index1, cut_index2, ss.npts


    # light curve
        
    # first compute the fluxes -- since we are going to normalize,
    # don't worry about the pi and sigma (Stefan-Boltzmann constant)

    # here we are assuming that R2 < R1
    R1 = float(R1)
    R2 = float(R2)

    f1 = float(R1**2 * T1**4)
    f2 = float(R2**2 * T2**4)

    f_normal = f1 + f2
    f_star2_transit = (R1**2 - R2**2) * T1**4 + R2**2 * T2**4
    f_star2_blocked = f1

    # relative fluxes -- normalized to normal
    f_star2_transit = f_star2_transit / f_normal
    f_star2_blocked = f_star2_blocked / f_normal

    print "f_star2_transit = ", f_star2_transit
    print "f_star2_blocked = ", f_star2_blocked

    # determine the times of the eclipses / transits
    # t_a = star 2 begins to pass in front of star 1
    # t_b = star 2 fully in front of star 1
    # t_c = star 2 begins to finish its transit of star 1
    # t_d = star 2 fully finished with transit
    # t_e = star 2 begins to go behind star 1
    # t_f = star 2 fully behind star 1
    # t_g = star 2 begins to emerge from behind star 1
    # t_h = star 2 fully emerged from behind star 1
    t_a = t_b = t_c = t_d = t_e = t_f = t_g = t_h = -1.0

    n = 0
    while (n < ss.npts):

        if (y_rel[n] <= 0):
            
            # star 2 in front of star 1
            if (x_rel[n] + R2*rad_scal > -R1*rad_scal and t_a == -1):
                t_a = ss.t[n]

            if (x_rel[n] - R2*rad_scal >= -R1*rad_scal and t_b == -1):
                t_b = ss.t[n]

            if (x_rel[n] + R2*rad_scal > R1*rad_scal and t_c == -1):
                t_c = ss.t[n]

            if (x_rel[n] - R2*rad_scal >= R1*rad_scal and t_d == -1):
                t_d = ss.t[n]

        else:

            # star 2 behind star 1
            if (x_rel[n] - R2*rad_scal < R1*rad_scal and t_e == -1):
                t_e = ss.t[n]

            if (x_rel[n] + R2*rad_scal <= R1*rad_scal and t_f == -1):
                t_f = ss.t[n]

            if (x_rel[n] - R2*rad_scal < -R1*rad_scal and t_g == -1):
                t_g = ss.t[n]

            if (x_rel[n] + R2*rad_scal <= -R1*rad_scal and t_h == -1):
                t_h = ss.t[n]
            

        n += 1


    # make an array of the flux vs. time -- this is the light curve
    f_system = numpy.zeros(ss.npts, numpy.float64)

    n = 0
    while (n < ss.npts):

        f_system[n] = 1.0
        
        # star 2 passing in front of star 1

        if (ss.t[n] >= t_a and ss.t[n] < t_b):
            # linearly interpolate between f = 1 at t = t_a and
            # f = f_star2_transit at t = t_b
            slope = (f_star2_transit - 1.0)/(t_b - t_a)
            f_system[n] = slope*(ss.t[n] - t_b) + f_star2_transit

        elif (ss.t[n] >= t_b and ss.t[n] < t_c):
            f_system[n] = f_star2_transit

        elif (ss.t[n] >= t_c and ss.t[n] < t_d):
            # linearly interpolate between f = f_star2_transit at 
            # t = t_c and f = 1 at t = t_d
            slope = (1.0 - f_star2_transit)/(t_d - t_c)
            f_system[n] = slope*(ss.t[n] - t_d) + 1.0

        # star 2 passing behind star 1

        elif (ss.t[n] >= t_e and ss.t[n] < t_f):
            # linearly interpolate between f = 1 at t = t_e and
            # f = f_star2_blocked at t = t_f
            slope = (f_star2_blocked - 1.0)/(t_f - t_e)
            f_system[n] = slope*(ss.t[n] - t_f) + f_star2_blocked

        elif (ss.t[n] >= t_f and ss.t[n] < t_g):
            f_system[n] = f_star2_blocked

        elif (ss.t[n] >= t_g and ss.t[n] < t_h):
            # linearly interpolate between f = f_star2_blocked 
            # at t = t_g and f = 1 at t = t_h
            slope = (1.0 - f_star2_blocked)/(t_h - t_g)
            f_system[n] = slope*(ss.t[n] - t_h) + 1.0

        n += 1


    # ================================================================
    # plotting
    # ================================================================


    iframe = 0

    n = 0
    while (n < ss.npts):

        pylab.clf()

        pylab.subplots_adjust(left=0.1,right=0.9,bottom=0.1,top=0.9)

        pylab.subplot(211)

        a = pylab.gca()
        a.set_aspect("equal", "datalim")
        pylab.axis("off")

        if (y_rel[n] < 0.0):
            # star 2 is in front of star 1

            # plot star 1's orbit position -- set to be the origin
            xc1, yc1 = circle(0, 0, R1*rad_scal)
            pylab.fill(xc1, yc1, 'r', ec="none")

            # plot star 2's relative orbit and position
            xc2, yc2 = circle(x_rel[n], y_rel[n], R2*rad_scal)
            pylab.fill(xc2, yc2, 'b', ec="none")

        else:
            # star 1 is in front of star 2

            # plot star 2's relative orbit and position
            xc2, yc2 = circle(x_rel[n], y_rel[n], R2*rad_scal)
            pylab.fill(xc2, yc2, 'b', ec="none")

            # plot star 1's orbit position -- set to be the origin
            xc1, yc1 = circle(0, 0, R1*rad_scal)
            pylab.fill(xc1, yc1, 'r', ec="none")



        # plot the orbit -- in two segment
        pylab.plot(x_rel[0:cut_index1],y_rel[0:cut_index1], 
                   color="0.75", linestyle="--",alpha=0.5)
        pylab.plot(x_rel[cut_index2:ss.npts],y_rel[cut_index2:ss.npts], 
                   color="0.75", linestyle="--",alpha=0.5)


        pylab.text(-1.5*a_star2, 1.1*a_star2, "star 1: T = %5.0f K" % (T1), 
                    color="r", fontsize=10)
        pylab.text(-1.5*a_star2, 0.95*a_star2, "star 2: T = %5.0f K" % (T2), 
                    color="b", fontsize=10)
        pylab.text(-1.5*a_star2, 0.8*a_star2, 
                    "R (star 1) / R (star 2) = %4.1f" % (R1/R2), 
                    color="k", fontsize=10)

        pylab.axis([-1.5*a_star2,1.5*a_star2,-1.25*a_star2,1.25*a_star2])


        pylab.subplot(212)

        pylab.plot(ss.t[0:ss.npts]/tmax, f_system, "k")
        pylab.scatter([ss.t[n]/tmax], [f_system[n]], color="k")

        pylab.xlim(0.0,1.0)
        pylab.ylim(0.0,1.1)
    
        pylab.xlabel("t/P")
        pylab.ylabel("relative flux")


        f = pylab.gcf()
        f.set_size_inches(5.0,6.5)
  
        outfile = "eclipsingbinar_%04d.png" % iframe
        pylab.savefig(outfile)


        iframe += 1
        n += 1




def circle(xc, yc, r):
    
    theta = numpy.arange(361)*2.0*math.pi/360

    x = r*numpy.cos(theta) + xc
    y = r*numpy.sin(theta) + yc

    return x, y



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


    
        
