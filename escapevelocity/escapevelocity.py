import math
import numpy
import pylab

# Shoot a projective horizontally some distance above Earth at various
# speeds and watch the resulting orbit
#
# Important speeds are the circular velocity (v_circ) and the escape
# velocity (v_escp)
#
# We consider:
#             v < v_circ   (crashes into Earth)
#             v = v_circ   (perfectly circular orbit)
#    v_circ < v < v_escp   (elliptical orbit)
#             v > v_escp   (escapes Earth's gravity)


# M. Zingale (2008-09-05)

# we work in CGS units
G = 6.67e-8
M_E = 5.9742e27
R_E = 6.378e8


# a simple class to serve as a container for the orbital information
class trajectory:
    
    def __init__ (self):

        self.npts = -1
        self.maxpoints = 2000
        self.x = numpy.zeros(self.maxpoints)
        self.y = numpy.zeros(self.maxpoints)
        self.t = numpy.zeros(self.maxpoints)




def orbit():

    # set the height from which we launch things
    hinit = 1.5*R_E
    max_rad = 25*R_E

    # the circular orbit will be our baseline
    orbit_circ = trajectory()

    vinit_circ = math.sqrt(G*M_E/hinit)
    tmax_circ = 2*math.pi*hinit/vinit_circ

    # we want dt to be constant for all of our orbits so we can compare
    # the velocities in the animation
    dt = tmax_circ/180.0        
    
    integrate_projectile(orbit_circ,vinit_circ,hinit,dt,max_rad)
    print "circ: ", orbit_circ.npts

    # v < v_c
    orbit_slow = trajectory()
    vinit_slow = 0.8*vinit_circ

    integrate_projectile(orbit_slow,vinit_slow,hinit,dt,max_rad)
    print "slow: ", orbit_slow.npts

    # v_c < v < v_e
    orbit_fast= trajectory()
    vinit_fast = 1.2*vinit_circ

    integrate_projectile(orbit_fast,vinit_fast,hinit,dt,max_rad)
    print "fast: ", orbit_fast.npts


    # v > v_e
    orbit_escp = trajectory()
    vinit_escp = 1.5*vinit_circ  # v_escp = sqrt(2) * v_circ

    integrate_projectile(orbit_escp,vinit_escp,hinit,dt,max_rad)
    print "escp: ", orbit_escp.npts


    # ================================================================
    # plotting
    # ================================================================

    # turn on interactive mode 
    #pylab.ion()
    img = pylab.imread("earth.png")


    # plot the orbits one by one
    iframe = 0

    # v < v_c
    n = 0
    while (n < orbit_slow.npts):

        pylab.clf()
        pylab.title(r"$v < v_\mathrm{circular}$", color="g", fontsize=20)

        # draw Earth -- we will use units that are in terms of Earth radii
        pylab.text(-0.95*4*R_E,-0.9*4*R_E, "Earth image credit:",
                    fontsize=7, color="0.50")
        pylab.text(-0.95*4*R_E,-0.95*4*R_E, "NASA/Apollo 17",
                    fontsize=7, color="0.50")
        pylab.imshow(img,extent = [-R_E,R_E,-R_E,R_E])

        # plot the current orbit
        pylab.plot(orbit_slow.x[0:n],orbit_slow.y[0:n], color="g")

        # print the time
        pylab.text(0.0,-0.95*4*R_E, "time = %6.4f hrs." % 
                   (orbit_slow.t[n]/3600.))

        pylab.axis([-4*R_E,4*R_E,-4*R_E,2*R_E])
        
        pylab.axis("off")

        f = pylab.gcf()
        f.set_size_inches(6.0,4.5)

        outfile = "escapevel_%04d.png" % iframe
        pylab.savefig(outfile)

        n += 1
        iframe += 1


    # v = v_c
    n = 0
    while (n < orbit_circ.npts):

        pylab.clf()
        pylab.title(r"$v = v_\mathrm{circular}$", color="r", fontsize=20)

        # draw Earth
        pylab.imshow(img,extent = [-R_E,R_E,-R_E,R_E])
        pylab.text(-0.95*4*R_E,-0.9*4*R_E, "Earth image credit:",
                    fontsize=7, color="0.50")
        pylab.text(-0.95*4*R_E,-0.95*4*R_E, "NASA/Apollo 17",
                    fontsize=7, color="0.50")

        # plot the previous orbit
        pylab.plot(orbit_slow.x[0:orbit_slow.npts-1],
                   orbit_slow.y[0:orbit_slow.npts-1], color="g",alpha="0.33")

        # plot the current orbit
        pylab.plot(orbit_circ.x[0:n],orbit_circ.y[0:n], color="r")

        # print the time
        pylab.text(0.0,-0.95*4*R_E, "time = %6.4f hrs." % 
                   (orbit_circ.t[n]/3600.))

        pylab.axis([-4*R_E,4*R_E,-4*R_E,2*R_E])
        
        pylab.axis("off")

        f = pylab.gcf()
        f.set_size_inches(6.0,4.5)

        outfile = "escapevel_%04d.png" % iframe
        pylab.savefig(outfile)

        n += 1
        iframe += 1


    # v_c < v < v_e
    n = 0
    while (n < orbit_fast.npts):

        pylab.clf()
        pylab.title(r"$v_\mathrm{circular} < v < v_\mathrm{escape}$", 
                    color="b", fontsize=20)

        # draw Earth
        pylab.text(-0.95*4*R_E,-0.9*4*R_E, "Earth image credit:",
                    fontsize=7, color="0.50")
        pylab.text(-0.95*4*R_E,-0.95*4*R_E, "NASA/Apollo 17",
                    fontsize=7, color="0.50")
        pylab.imshow(img,extent = [-R_E,R_E,-R_E,R_E])

        # plot the previous orbits
        pylab.plot(orbit_slow.x[0:orbit_slow.npts-1],
                   orbit_slow.y[0:orbit_slow.npts-1], color="g", alpha="0.33")
        pylab.plot(orbit_circ.x[0:orbit_circ.npts-1],
                   orbit_circ.y[0:orbit_circ.npts-1], color="r", alpha="0.33")

        # plot the current orbit
        pylab.plot(orbit_fast.x[0:n],orbit_fast.y[0:n], color="b")

        # print the time
        pylab.text(0.0,-0.95*4*R_E, "time = %6.4f hrs." % 
                   (orbit_fast.t[n]/3600.))

        pylab.axis([-4*R_E,4*R_E,-4*R_E,2*R_E])
        
        pylab.axis("off")

        f = pylab.gcf()
        f.set_size_inches(6.0,4.5)

        outfile = "escapevel_%04d.png" % iframe
        pylab.savefig(outfile)

        n += 1
        iframe += 1


    # v > e_c
    n = 0
    while (n < orbit_escp.npts):

        pylab.clf()
        pylab.title(r"$v > v_\mathrm{escape}$", color="k", fontsize=20)
        pylab.imshow(img,extent = [-R_E,R_E,-R_E,R_E])

        # plot the previous orbits
        pylab.plot(orbit_slow.x[0:orbit_slow.npts-1],
                   orbit_slow.y[0:orbit_slow.npts-1], color="g", alpha="0.33")
        pylab.plot(orbit_circ.x[0:orbit_circ.npts-1],
                   orbit_circ.y[0:orbit_circ.npts-1], color="r", alpha="0.33")
        pylab.plot(orbit_fast.x[0:orbit_fast.npts-1],
                   orbit_fast.y[0:orbit_fast.npts-1], color="b", alpha="0.33")

        
        # plot the current orbit
        pylab.plot(orbit_escp.x[0:n+1],orbit_escp.y[0:n+1],color="k")


        # for this case, we want to zoom out in a smooth fashion
        max_coord = max(max(math.fabs(orbit_escp.x[n+1]),
                            math.fabs(orbit_escp.y[n+1])),
                        4*R_E)

        pylab.axis([-max_coord,max_coord,
                    -max_coord,0.5*max_coord])
        
        pylab.axis("off")

        f = pylab.gcf()
        f.set_size_inches(6.0,4.5)

        # print the image credit
        pylab.text(-0.95*max_coord,-0.9*max_coord, "Earth image credit:",
                    fontsize=7, color="0.50")
        pylab.text(-0.95*max_coord,-0.95*max_coord, "NASA/Apollo 17",
                    fontsize=7, color="0.50")

        # print the time
        pylab.text(0.0,-0.95*max_coord, "time = %6.4f hrs." % 
                   (orbit_escp.t[n]/3600.))

        outfile = "escapevel_%04d.png" % iframe
        pylab.savefig(outfile)

        n += 1
        iframe += 1





def integrate_projectile(orbit,vinit,hinit,dt,max_rad):

    SMALL = 1.e-16

    # allocate storage for R-K intermediate results
    k1 = numpy.zeros(8, numpy.float64)
    k2 = numpy.zeros(8, numpy.float64)
    k3 = numpy.zeros(8, numpy.float64)
    k4 = numpy.zeros(8, numpy.float64)


    y = numpy.zeros(8, numpy.float64)
    f = numpy.zeros(8, numpy.float64)



    t = 0.0

    # initial conditions
    y[0] = 0       # initial x position
    y[1] = hinit   # initial height

    y[2] = vinit   # initial tangental velocity
    y[3] = 0.0     # initial vertical (radial) velocity

    orbit.x[0] = y[0]
    orbit.y[0] = y[1]
    orbit.t[0] = t

    r = math.sqrt(y[0]**2 + y[1]**2)


    numorbits = 0
    angle_old = math.pi/2   # we are initially vertical

    n = 1
    while (n < orbit.maxpoints and r > R_E and r <= max_rad and numorbits == 0):
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

        orbit.x[n] = y[0]
        orbit.y[n] = y[1]
        orbit.t[n] = t

        # compute the radius to decide if we collided with Earth
        r = math.sqrt(y[0]**2 + y[1]**2)

        # compute the angle wrt to vertical and compare to the 
        # previous one.  If we go from > pi/2 to <= pi/2 then we have
        # completed an orbit
        angle = math.atan2(y[1],(y[0] + SMALL))

        if (angle_old > math.pi/2 and angle <= math.pi/2):
            numorbits = 1

        angle_old = angle

        n += 1

    
    orbit.npts = n



def rhs(t,y):

    f = numpy.zeros(8, numpy.float64)

    # y[0] = x, y[1] = y, y[2] = v_x, y[3] = v_y
    f[0] = y[2]
    f[1] = y[3]
    r = numpy.sqrt(y[0]**2 + y[1]**2)
    f[2] = -G*M_E*y[0]/r**3
    f[3] = -G*M_E*y[1]/r**3

    return f
    

    
if __name__== "__main__":
    orbit()


    
        
