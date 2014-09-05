import math
import numpy as np
import matplotlib.pyplot as plt
import anim_solvers.earth_orbit as earth_orbit

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


def orbit():

    # set the height from which we launch things
    hinit = 1.5*R_E
    max_rad = 25*R_E

    # the circular orbit will be our baseline
    orbit_circ = earth_orbit.trajectory(GM=G*M_E, R_crash=R_E)

    vinit_circ = math.sqrt(G*M_E/hinit)
    tmax_circ = 2*math.pi*hinit/vinit_circ

    # we want dt to be constant for all of our orbits so we can compare
    # the velocities in the animation
    dt = tmax_circ/180.0        
    
    orbit_circ.integrate(vinit_circ, hinit, dt, max_rad)
    print "circ: ", orbit_circ.npts

    # v < v_c
    orbit_slow = earth_orbit.trajectory(GM=G*M_E, R_crash=R_E)
    vinit_slow = 0.8*vinit_circ

    orbit_slow.integrate(vinit_slow, hinit, dt, max_rad)
    print "slow: ", orbit_slow.npts

    # v_c < v < v_e
    orbit_fast= earth_orbit.trajectory(GM=G*M_E, R_crash=R_E)
    vinit_fast = 1.2*vinit_circ

    orbit_fast.integrate(vinit_fast, hinit, dt, max_rad)
    print "fast: ", orbit_fast.npts


    # v > v_e
    orbit_escp = earth_orbit.trajectory(GM=G*M_E, R_crash=R_E)
    vinit_escp = 1.5*vinit_circ  # v_escp = sqrt(2) * v_circ

    orbit_escp.integrate(vinit_escp, hinit, dt, max_rad)
    print "escp: ", orbit_escp.npts


    # ================================================================
    # plotting
    # ================================================================

    img = plt.imread("earth.png")


    # plot the orbits one by one
    iframe = 0

    # v < v_c
    for n in range(orbit_slow.npts):

        plt.clf()
        plt.title(r"$v < v_\mathrm{circular}$", color="g", fontsize=20)

        # draw Earth -- we will use units that are in terms of Earth radii
        plt.imshow(img,extent = [-R_E,R_E,-R_E,R_E])

        # plot the current orbit
        plt.plot(orbit_slow.x[0:n],orbit_slow.y[0:n], color="g")

        plt.axis([-4*R_E,4*R_E,-4*R_E,2*R_E])
        
        plt.axis("off")

        ax = plt.gca()
        ax.set_aspect("equal", "datalim")

        f = plt.gcf()
        f.set_size_inches(9.6,7.2)

        plt.text(0.05, 0.06, "Earth image credit:", transform=f.transFigure,
                 fontsize=7, color="0.50")
        plt.text(0.05, 0.04, "NASA/Apollo 17", transform=f.transFigure,
                    fontsize=7, color="0.50")

        # print the time
        plt.text(0.7,0.1, "time = %6.4f hrs." % (orbit_slow.t[n]/3600.),
                 transform=f.transFigure)


        plt.savefig("escapevel_%04d.png" % iframe)

        iframe += 1


    # v = v_c
    for n in range(orbit_circ.npts):

        plt.clf()
        plt.title(r"$v = v_\mathrm{circular}$", color="r", fontsize=20)

        # draw Earth
        plt.imshow(img,extent = [-R_E,R_E,-R_E,R_E])

        # plot the previous orbit
        plt.plot(orbit_slow.x[0:orbit_slow.npts-1],
                   orbit_slow.y[0:orbit_slow.npts-1], color="g",alpha=0.33)

        # plot the current orbit
        plt.plot(orbit_circ.x[0:n],orbit_circ.y[0:n], color="r")

        plt.axis([-4*R_E,4*R_E,-4*R_E,2*R_E])
        
        plt.axis("off")

        ax = plt.gca()
        ax.set_aspect("equal", "datalim")

        f = plt.gcf()
        f.set_size_inches(9.6,7.2)

        plt.text(0.05, 0.06, "Earth image credit:", transform=f.transFigure,
                 fontsize=7, color="0.50")
        plt.text(0.05, 0.04, "NASA/Apollo 17", transform=f.transFigure,
                    fontsize=7, color="0.50")

        # print the time 
        plt.text(0.7,0.1, "time = %6.4f hrs." % (orbit_circ.t[n]/3600.),
                 transform=f.transFigure)

        plt.savefig("escapevel_%04d.png" % iframe)

        iframe += 1


    # v_c < v < v_e
    for n in range(orbit_fast.npts):

        plt.clf()
        plt.title(r"$v_\mathrm{circular} < v < v_\mathrm{escape}$", 
                    color="b", fontsize=20)

        # draw Earth
        plt.imshow(img,extent = [-R_E,R_E,-R_E,R_E])

        # plot the previous orbits
        plt.plot(orbit_slow.x[0:orbit_slow.npts-1],
                   orbit_slow.y[0:orbit_slow.npts-1], color="g", alpha=0.33)
        plt.plot(orbit_circ.x[0:orbit_circ.npts-1],
                   orbit_circ.y[0:orbit_circ.npts-1], color="r", alpha=0.33)

        # plot the current orbit
        plt.plot(orbit_fast.x[0:n],orbit_fast.y[0:n], color="b")

        plt.axis([-4*R_E,4*R_E,-4*R_E,2*R_E])
        
        plt.axis("off")

        ax = plt.gca()
        ax.set_aspect("equal", "datalim")

        f = plt.gcf()
        f.set_size_inches(9.6,7.2)

        plt.text(0.05, 0.06, "Earth image credit:", transform=f.transFigure,
                 fontsize=7, color="0.50")
        plt.text(0.05, 0.04, "NASA/Apollo 17", transform=f.transFigure,
                    fontsize=7, color="0.50")

        # print the time
        plt.text(0.7,0.1, "time = %6.4f hrs." % (orbit_fast.t[n]/3600.),
                 transform=f.transFigure)

        plt.savefig("escapevel_%04d.png" % iframe)

        iframe += 1


    # v > e_c
    for n in range(orbit_escp.npts):

        plt.clf()
        plt.title(r"$v > v_\mathrm{escape}$", color="k", fontsize=20)
        plt.imshow(img,extent = [-R_E,R_E,-R_E,R_E])

        # plot the previous orbits
        plt.plot(orbit_slow.x[0:orbit_slow.npts-1],
                   orbit_slow.y[0:orbit_slow.npts-1], color="g", alpha=0.33)
        plt.plot(orbit_circ.x[0:orbit_circ.npts-1],
                   orbit_circ.y[0:orbit_circ.npts-1], color="r", alpha=0.33)
        plt.plot(orbit_fast.x[0:orbit_fast.npts-1],
                   orbit_fast.y[0:orbit_fast.npts-1], color="b", alpha=0.33)

        
        # plot the current orbit
        plt.plot(orbit_escp.x[0:n+1],orbit_escp.y[0:n+1],color="k")


        # for this case, we want to zoom out in a smooth fashion
        max_coord = max(max(math.fabs(orbit_escp.x[n+1]),
                            math.fabs(orbit_escp.y[n+1])),
                        4*R_E)

        plt.axis([-max_coord,max_coord,
                    -max_coord,0.5*max_coord])
        
        plt.axis("off")

        ax = plt.gca()
        ax.set_aspect("equal", "datalim")

        f = plt.gcf()
        f.set_size_inches(9.6,7.2)

        plt.text(0.05, 0.06, "Earth image credit:", transform=f.transFigure,
                 fontsize=7, color="0.50")
        plt.text(0.05, 0.04, "NASA/Apollo 17", transform=f.transFigure,
                    fontsize=7, color="0.50")

        # print the time
        plt.text(0.7,0.1, "time = %6.4f hrs." % (orbit_escp.t[n]/3600.),
                 transform=f.transFigure)

        plt.savefig("escapevel_%04d.png" % iframe)

        iframe += 1


    
if __name__== "__main__":
    orbit()


    
        
