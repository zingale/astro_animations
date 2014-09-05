import math
import numpy as np
import matplotlib.pyplot as plt
import anim_solvers.earth_orbit as earth_orbit

# Shoot a projective horizontally some distance above Earth at various
# speeds and watch the resulting orbit
#
# Give it a boost a perihelion and observe the new orbit.

# M. Zingale (2008-09-05)

# we work in CGS units
G = 6.67e-8
M_E = 5.9742e27
R_E = 6.378e8
SMALL = 1.e-12

max_coord = 3*R_E

def orbit():

    # set the height from which we launch things
    hinit = 1.5*R_E
    max_rad = 25*R_E


    # the circular orbit will be our baseline
    vinit_circ = math.sqrt(G*M_E/hinit)
    tmax_circ = 2*math.pi*hinit/vinit_circ

    # we want dt to be constant for all of our orbits so we can compare
    # the velocities in the animation
    dt = tmax_circ/240.0        
    

    # orbit 1: v = v_c
    orbit1 = earth_orbit.trajectory(GM=G*M_E, R_crash=R_E)
    vinit1 = vinit_circ

    orbit1.integrate(vinit1, hinit, dt, max_rad)
    print "orbit 1: ", orbit1.npts

    # orbit 2: v = 1.1*v_c
    orbit2 = earth_orbit.trajectory(GM=G*M_E, R_crash=R_E)
    vinit2 = 1.1*vinit_circ

    orbit2.integrate(vinit2, hinit, dt, max_rad)
    print "orbit 2: ", orbit2.npts


    # orbit 3: v = (1.1)**2*v_c
    orbit3 = earth_orbit.trajectory(GM=G*M_E, R_crash=R_E)
    vinit3 = 1.1*1.1*vinit_circ

    orbit3.integrate(vinit3, hinit, dt, max_rad)
    print "orbit 3: ", orbit3.npts

    print np.min(orbit3.y)/R_E
    print np.max(orbit3.y)/R_E

    # ================================================================
    # plotting
    # ================================================================

    img = plt.imread("earth.png")

    # plot the orbits one by one
    iframe = 0

    # orbit 1 (circular)

    for n in range(orbit1.npts):

        plt.clf()
 
        plt.title(r"$\mathrm{circular\ orbit:}\ v_\mathrm{peri} = v_\mathrm{circular}$", color="g", fontsize=16)

        # draw Earth -- we will use units that are in terms of Earth radii
        plt.imshow(img,extent = [-R_E,R_E,-R_E,R_E])

        # plot the current orbit
        plt.plot(orbit1.x[0:n+1],orbit1.y[0:n+1], color="g")

        # draw the spaceship
        plt.scatter([orbit1.x[n]],[orbit1.y[n]], color="b", marker="o")

        angle = math.atan2(orbit1.y[n],(orbit1.x[n] + SMALL))*180/math.pi

        # draw flames for the boost at the start of next orbit
        if (angle < 102.5 and angle >= 89 and n > 0.75*orbit1.npts):
            plt.scatter([orbit1.x[n]-0.025*max_coord],[orbit1.y[n]], color="r", marker="<")


        plt.axis("off")

        ax = plt.gca()
        ax.set_aspect("equal")

        f = plt.gcf()
        f.set_size_inches(9.6,7.2)

        plt.xlim(-3*R_E, 3*R_E)
        plt.ylim(-5*R_E, 3*R_E)
                

        # print the image credit
        plt.text(0.05, 0.06, "Earth image credit:", transform=f.transFigure,
                 fontsize=7, color="0.50")
        plt.text(0.05, 0.04, "NASA/Apollo 17", transform=f.transFigure,
                    fontsize=7, color="0.50")

        # print the time
        plt.text(0.7,0.1, "time = %6.4f hrs." % (orbit1.t[n]/3600.),
                 transform=f.transFigure)

        plt.savefig("changeorbit_%04d.png" % iframe)

        iframe += 1


    # orbit 2 (first boost)

    for n in range(orbit2.npts):

        plt.clf()

        plt.title(r"$\mathrm{first\ boost:}\ v_\mathrm{peri} = 1.1 \times\ v_\mathrm{circular}$", color="g", fontsize=16)

        # draw Earth
        plt.imshow(img,extent = [-R_E,R_E,-R_E,R_E])

        # plot the previous orbit
        plt.plot(orbit1.x[0:orbit1.npts],
                   orbit1.y[0:orbit1.npts], color="k",alpha=0.33)

        # plot the current orbit
        plt.plot(orbit2.x[0:n+1],orbit2.y[0:n+1], color="g")

        # draw the spaceship
        plt.scatter([orbit2.x[n]],[orbit2.y[n]], color="b", marker="o")

        angle = math.atan2(orbit2.y[n],(orbit2.x[n] + SMALL))*180/math.pi

        # draw flames for the boost at start of this orbit
        if (angle <= 90 and angle > 78.5 and n < 0.25*orbit2.npts):
            plt.scatter([orbit2.x[n]-0.025*max_coord],[orbit2.y[n]], color="r", marker="<")
        
        # draw flames for the boost at the start of next orbit
        if (angle < 102.5 and angle >= 89 and n > 0.75*orbit2.npts):
            plt.scatter([orbit2.x[n]-0.025*max_coord],[orbit2.y[n]], color="r", marker="<")
        
        plt.axis([-3*R_E, 3*R_E, -5*R_E, 3*R_E])

        plt.axis("off")

        ax = plt.gca()
        ax.set_aspect("equal")

        f = plt.gcf()
        f.set_size_inches(9.6,7.2)

        plt.text(0.05, 0.06, "Earth image credit:", transform=f.transFigure,
                 fontsize=7, color="0.50")
        plt.text(0.05, 0.04, "NASA/Apollo 17", transform=f.transFigure,
                    fontsize=7, color="0.50")

        # print the time
        plt.text(0.7,0.1, "time = %6.4f hrs." % (orbit2.t[n]/3600.),
                 transform=f.transFigure)

        plt.savefig("changeorbit_%04d.png" % iframe)

        iframe += 1


    # orbit 3 (second boost)
    for n in range(orbit3.npts):

        plt.clf()

        plt.title(r"$\mathrm{second\ boost:}\ v_\mathrm{peri} = (1.1)^2 \times\ v_\mathrm{circular}$", color="g", fontsize=16)

        # draw Earth
        plt.imshow(img,extent = [-R_E,R_E,-R_E,R_E])


        # plot the previous orbit
        plt.plot(orbit1.x[0:orbit1.npts],
                   orbit1.y[0:orbit1.npts], color="k",alpha=0.33)

        plt.plot(orbit2.x[0:orbit2.npts],
                   orbit2.y[0:orbit2.npts], color="k",alpha=0.33)


        # plot the current orbit
        plt.plot(orbit3.x[0:n+1],orbit3.y[0:n+1], color="g")

        # draw the spaceship
        plt.scatter([orbit3.x[n]],[orbit3.y[n]], color="b", marker="o")


        angle = math.atan2(orbit3.y[n],(orbit3.x[n] + SMALL))*180/math.pi

        # draw flames for the boost at the start of this orbit
        if (angle <= 90 and angle > 78.5 and n < 0.25*orbit3.npts):
            plt.scatter([orbit3.x[n]-0.025*max_coord],[orbit3.y[n]], color="r", marker="<")

        # draw flames for the boost at the start of next orbit
        if (angle >= -90 and angle < -78.5 and n > 0.75*orbit3.npts):
            plt.scatter([orbit3.x[n]-0.025*max_coord],[orbit3.y[n]], color="r", marker="<")


        plt.axis([-3*R_E, 3*R_E, -5*R_E, 3*R_E])
          
        plt.axis("off")

        ax = plt.gca()
        ax.set_aspect("equal")

        f = plt.gcf()
        f.set_size_inches(9.6,7.2)

        plt.text(0.05, 0.06, "Earth image credit:", transform=f.transFigure,
                 fontsize=7, color="0.50")
        plt.text(0.05, 0.04, "NASA/Apollo 17", transform=f.transFigure,
                    fontsize=7, color="0.50")

        # print the time
        plt.text(0.7,0.1, "time = %6.4f hrs." % (orbit3.t[n]/3600.),
                 transform=f.transFigure)

        plt.savefig("changeorbit_%04d.png" % iframe)

        iframe += 1

    
if __name__== "__main__":
    orbit()


    
        
