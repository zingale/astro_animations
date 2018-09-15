import numpy as np
import matplotlib.pyplot as plt

# show the ellipse geometry

# M. Zingale (2009-02-13)


def ellipse():

    # theta ranges from 0 to 2pi
    npts = 250
    theta = np.arange(npts)*2.0*np.pi/(npts-1)

    e = 0.5
    a = 1.0
    b = a*np.sqrt(1-e*e)

    n = npts//5

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

    # draw our current point
    plt.scatter([x[n]],[y[n]], color="k", s=25)

    # second foci
    plt.scatter([-2.0*a*e],[0], color="g", marker="x", s=200)

    # primary foci
    plt.scatter([0], [0], color="g", marker="x", s=200)

    # center
    plt.scatter([-a*e], [0], color='k', marker="x", s=200)

    # draw the semi-major axis
    plt.plot([-a*e,a*(1-e)],[0,0], color='r')
    plt.text(0.5*a - a*e, -0.15*b, "a", color='r')

    # draw the semi-minor axis
    plt.plot([-a*e,-a*e], [0,-b], color='k')
    plt.text(-a*e-0.15*a*e, -0.5*b, "b", color='k')

    plt.plot([-a*e, -2*a*e], [0,0], color="0.5")
    plt.text(-1.5*a*e, -0.15*b, "ae", color="0.5")

    # draw lines connecting the foci to the current point
    plt.plot([0,x[n]],[0,y[n]], color="g")
    plt.text(0, 0.3*b, "r", color="g")

    plt.plot([-2.0*a*e,x[n]],[0,y[n]], color="g")
    plt.text(-a*e-0.1*a, 0.3*b, r"r$^\prime$", color="g")
    # indicate the angle
    plt.text(0.075*a, 0.05*b, r"$\theta$", color="k")

    plt.axis([-1.75,.75,-1.25,1.25])

    f = plt.gcf()
    f.set_size_inches(6.0,6.0)

    outfile = "ellipse_geom.png"
    plt.savefig(outfile)


if __name__== "__main__":
    ellipse()
