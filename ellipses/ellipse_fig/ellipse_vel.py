import numpy as np
import matplotlib.pyplot as plt

# show the velocity vectors at perihelion and aphelion for ellipse geometry

# M. Zingale (2010-09-08)


def ellipse():

    # theta ranges from 0 to 2pi
    npts = 250
    theta = np.arange(npts)*2.0*np.pi/(npts-1)

    e = 0.5
    a = 1.0
    b = a*np.sqrt(1-e*e)

    n = npts/5

    r = a*(1.0 - e*e)/(1.0 + e*np.cos(theta))

    x = r*np.cos(theta)
    y = r*np.sin(theta)

    # plotting
    plt.clf()

    ax = plt.gca()
    ax.set_aspect("equal", "datalim")
    plt.axis("off")

    # Draw the ellipse
    plt.plot(x,y,color="b")

    # primary foci
    plt.scatter([0], [0], color="k", marker="x", s=200)

    # draw arrow to perihelion
    eps = 0.033*a
    plt.plot([0,a*(1-e)],[0,0], color='r')
    plt.plot([a*(1-e),a*(1-e)-eps],[0,eps], color='r')
    plt.plot([a*(1-e),a*(1-e)-eps],[0,-eps], color='r')
    plt.text(0.5*a*(1-e), -0.15*b, "r", color='r')


    # draw velocity at perihelion
    vel = 0.66*a
    plt.plot([a*(1-e),a*(1-e)], [0, vel], color='g')
    plt.plot([a*(1-e),a*(1-e)-eps], [vel, vel-eps], color='g')
    plt.plot([a*(1-e),a*(1-e)+eps], [vel, vel-eps], color='g')
    plt.text(a*(1-e)+0.05*b,0.5*vel, "v", color='g')

    plt.axis([-1.75,.75,-1.25,1.25])


    f = plt.gcf()
    f.set_size_inches(6.0,6.0)

    outfile = "ellipse_vel.png"
    plt.savefig(outfile)


if __name__== "__main__":
    ellipse()
