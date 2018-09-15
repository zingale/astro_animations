import numpy as np
import matplotlib.pyplot as plt

# show the initial conditions for integrating an orbit

def ellipse():

    # theta ranges from 0 to 2pi
    npts = 250
    theta = np.arange(npts)*2.0*np.pi/(npts-1)

    e = 0.5
    a = 1.0
    b = a*np.sqrt(1-e*e)

    n = npts/5

    r = a*(1.0 - e*e)/(1.0 + e*np.cos(theta))

    # perihelion at +y
    x = r*np.cos(theta + np.pi/2)
    y = r*np.sin(theta + np.pi/2)

    # plotting
    plt.clf()

    ax = plt.gca()
    ax.set_aspect("equal", "datalim")
    plt.axis("off")

    # draw the ellipse
    plt.plot(x, y, color="C0")

    # primary foci
    plt.scatter([0], [0], color="k", marker=(20,1), s=250, zorder=-100)
    plt.scatter([0], [0], color="y", marker=(20,1), s=200)

    xp = 0
    yp = a*(1-e)

    # center
    plt.scatter([0], [-a*e], color='k', marker="x", s=50)

    # semi-major axis
    plt.plot([0, 0], [-a*e, -a*(1+e)], color="0.5")
    plt.text(-0.1*b, -0.5*(a*e + a*(1+e)), "a", color="0.5")

    # draw arrow to perihelion
    eps = 0.033*a
    plt.plot([0, 0], [0, yp], color='C1')
    plt.plot([0, eps], [a*(1-e), a*(1-e)-eps], color='C1')
    plt.plot([0, -eps], [a*(1-e), a*(1-e)-eps], color='C1')
    plt.text(-0.1*b, 0.5*a*(1-e), "r", color='C1')


    # draw velocity at perihelion
    vel = 0.5*a
    plt.plot([0, -vel], [a*(1-e), a*(1-e)], color='C2')
    plt.plot([-vel, -vel+eps], [a*(1-e), a*(1-e)-eps], color='C2')
    plt.plot([-vel, -vel+eps], [a*(1-e), a*(1-e)+eps], color='C2')
    plt.text(-0.5*vel, a*(1-e)+0.05*b, "v", color='C2')

    plt.axis([-1.25, 1.25, -1.75, 0.75])

    f = plt.gcf()
    f.set_size_inches(6.0,6.0)

    outfile = "ellipse_initial_cond.png"
    plt.tight_layout()
    plt.savefig(outfile, dpi=150, bbox_inches="tight")


if __name__== "__main__":
    ellipse()
