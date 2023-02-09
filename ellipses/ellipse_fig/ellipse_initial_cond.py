import numpy as np
import matplotlib.pyplot as plt

# show the initial conditions for integrating an orbit

def ellipse():

    # theta ranges from 0 to 2pi
    npts = 360
    theta = np.arange(npts)*2.0*np.pi/(npts-1)

    # semi-major axis, eccentricity, and semi-minor axis
    e = 0.4
    a = 1.0
    b = a * np.sqrt(1 - e * e)

    # equation of an ellipse, where r is the distance
    # from a foci

    r = a * (1.0 - e * e) / (1.0 + e * np.cos(theta))

    # perihelion at +y
    x = r*np.cos(theta + np.pi/2)
    y = r*np.sin(theta + np.pi/2)

    fig, ax = plt.subplots()

    ax.set_aspect("equal", "datalim")
    ax.axis("off")

    # draw the ellipse
    ax.plot(x, y, color="C0")

    # primary foci
    ax.scatter([0], [0], color="k", marker=(20,1), s=250, zorder=-100)
    ax.scatter([0], [0], color="y", marker=(20,1), s=200)

    xp = 0
    yp = a * (1 - e)

    # center
    ax.scatter([0], [-a*e], color='k', marker="x", s=50)

    # semi-major axis
    ax.plot([0, 0], [-a*e, -a*(1+e)], color="0.5")
    ax.text(-0.1*b, -0.5*(a*e + a*(1+e)), "a", color="0.5")

    # draw arrow to perihelion
    eps = 0.033*a
    ax.plot([0, 0], [xp, yp], color='C1')
    ax.plot([0, eps], [a*(1-e), a*(1-e)-eps], color='C1')
    ax.plot([0, -eps], [a*(1-e), a*(1-e)-eps], color='C1')
    ax.text(-0.1*b, 0.5*a*(1-e), "r", color='C1')


    # draw velocity at perihelion
    vel = 0.5*a
    ax.plot([0, -vel], [a*(1-e), a*(1-e)], color='C2')
    ax.plot([-vel, -vel+eps], [a*(1-e), a*(1-e)-eps], color='C2')
    ax.plot([-vel, -vel+eps], [a*(1-e), a*(1-e)+eps], color='C2')
    ax.text(-0.5*vel, a*(1-e)+0.05*b, "v", color='C2')

    ax.axis([-1.25 * a, 1.25*a, -1.25 * a * (1 + e), 1.25 * a * (1 - e)])

    fig.set_size_inches(5.0, 5.0)

    outfile = f"ellipse_initial_cond_e{e}.png"
    fig.tight_layout()
    fig.savefig(outfile, dpi=200, bbox_inches="tight")


if __name__== "__main__":
    ellipse()
