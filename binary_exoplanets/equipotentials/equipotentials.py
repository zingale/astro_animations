# draw lines of equipotentials for a rotating binary system.
#
# The techniques here come from "Astrophysics with a PC: An
# Introduction to Computational Astrophysics" by Paul Hellings

import math
import numpy as np
import matplotlib.pyplot as plt

def equipotentials(mu):
    # here mu is the mass parameter.  mu = 1/2 is equal mass
    # this is related to the q parameter (M2 / M1; with q <= 1).  

    # note, if you make mu significant at many digits, you need to
    # have many points to resolve features
    npts = 1024

    xmin = -2.0
    xmax = 2.0
    ymin = -2.0
    ymax = 2.0

    EPS = 1.e-12

    NITER = 100

    x = np.linspace(xmin, xmax, npts, dtype=np.float64)
    y = np.linspace(ymin, ymax, npts, dtype=np.float64)

    X, Y = np.meshgrid(x, y)

    # this is the dimensionless potential, in the z=0 plane
    V = Vf(mu, X, Y)


    # find L1
    x0 = 0.0

    for n in range(NITER):
    
        dVX = dVXdx(-1.0, 1.0, mu, x0)
        d2VX = d2VXdx2(-1.0, 1.0, mu, x0)

        x1 = x0 - dVX/d2VX
        err = abs(x1 - x0)/abs(x0 + EPS)

        if err < EPS: break

        x0 = x1

    x_L1 = x1
    y_L1 = 0.0
    V_L1 = Vf(mu, x_L1, y_L1)


    # find L2
    x0 = -1.0

    for n in range(NITER):
    
        dVX = dVXdx(-1.0, -1.0, mu, x0)
        d2VX = d2VXdx2(-1.0, -1.0, mu, x0)

        x1 = x0 - dVX/d2VX
        err = abs(x1 - x0)/abs(x0 + EPS)

        if err < EPS: break

        x0 = x1


    x_L2 = x1
    y_L2 = 0.0
    V_L2 = Vf(mu, x_L2, y_L2)


    # find L3
    x0 = 1.0

    for n in range(NITER):

        dVX = dVXdx(1.0, 1.0, mu, x0)
        d2VX = d2VXdx2(1.0, 1.0, mu, x0)

        x1 = x0 - dVX/d2VX
        err = abs(x1 - x0)/abs(x0 + EPS)

        if err < EPS: break

        x0 = x1


    x_L3 = x1
    y_L3 = 0.0
    V_L3 = Vf(mu, x_L3, y_L3)

    # L4 and L5
    x_L4 = mu - 0.5
    y_L4 = math.sin(math.pi/3.0)

    x_L5 = mu - 0.5
    y_L5 = -math.sin(math.pi/3.0)

    
    # do imshow
    plt.clf()

    print mu,  V.min(), V.max()
    plt.imshow(np.log10(V), origin="lower", cmap="Accent",
               extent=[xmin, xmax, ymin, ymax])

    # draw contours -- these values seem reasonable for a range of mu's
    Vmin = 1.5 
    Vmax = 1000.0  # np.max(V)
    nC = 25
    
    C = np.logspace(math.log10(Vmin), math.log10(Vmax), nC)

    plt.contour(x, y, V, C, colors="b")

    # special contours right through the lagrange points
    plt.contour(x, y, V, [V_L1], colors="b")
    plt.contour(x, y, V, [V_L2], colors="b")
    plt.contour(x, y, V, [V_L3], colors="b")
 
    
    # mark the Lagrange points and write the names
    xeps = 0.025

    plt.scatter([x_L1], [y_L1], marker="x", color="r", s=50)
    plt.text(x_L1+xeps, y_L1+xeps, "L1", color="r")

    plt.scatter([x_L2], [y_L2], marker="x", color="r", s=50)
    plt.text(x_L2+xeps, y_L2+xeps, "L2", color="r")

    plt.scatter([x_L3], [y_L3], marker="x", color="r", s=50)
    plt.text(x_L3+xeps, y_L3+xeps, "L3", color="r")

    plt.scatter([x_L4], [y_L4], marker="x", color="r", s=50)
    plt.text(x_L4+xeps, y_L4+xeps, "L4", color="r")

    plt.scatter([x_L5], [y_L5], marker="x", color="r", s=50)
    plt.text(x_L5+xeps, y_L5+xeps, "L5", color="r")
       
    plt.axis([xmin,xmax,ymin,ymax])

    plt.title(r"Equipotentials, mu = {:5.3f}".format(mu) ,fontsize=12)
    plt.xlabel("x/(a + b)")
    plt.ylabel("y/(a + b)")

    f = plt.gcf()
    f.set_size_inches(10.8,10.8)

    plt.tight_layout()

    plt.savefig("equipotentials_mu_{:5.3f}.png".format(mu))


def Vf(mu, x, y):
    V = (1.0 - mu)/np.sqrt( (x - mu)**2 + y**2 ) + \
        mu/np.sqrt( (x + 1.0 - mu)**2 + y**2 ) + \
        0.5*(x**2 + y**2)    
    return V


def dVXdx(h1, h2, mu, x):
    dVX = -h1*(1.0 - mu)/(x - mu)**2 - h2*mu/(x + 1.0 - mu)**2 + x
    return dVX


def d2VXdx2(h1, h2, mu, x):
    d2VX = 2.0*h1*(1.0 - mu)/(x - mu)**3 + 2.0*h2*mu/(x + 1.0 - mu)**3 + 1.0
    return d2VX



if __name__== "__main__":

    for mu in np.linspace(0.005, 0.5, 200):
        equipotentials(mu)


