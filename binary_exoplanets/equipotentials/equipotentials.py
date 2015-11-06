# draw lines of equipotentials for a rotating binary system.
#
# The techniques here come from "Astrophysics with a PC: An
# Introduction to Computational Astrophysics" by Paul Hellings

import math
import numpy as np
import matplotlib.pyplot as plt



class Equipotentials(object):
    """the general setup is two stars on the x-axis.  M_1 is the more
    massive and is at x = +a.  M_2 is the less massive and is at
    x = -b.
    
    center of mass tells us that: M_1 a = M_2 b
    
    we work in mass units of (M_1 + M_2), and distance units of 
    (a + b)
    
    here mu is the mass parameter, which we can think of as 

    mu = M_2 / (M_1 + M_2).  

    mu = 1/2 is equal mass, but in general, with M_1 > M_2, mu
    will be less than 1.  The normalized mass of the primary (more
    massive star) is then 1 - mu
    
    In these reduced coordinates, the distance a becomes 
    
    a / (a + b) = mu.  
    
    We put the center of mass at the origin, so this is at x =
    mu. (negative x axis) likewise star 2 is at 
    
    b / (a + b) = (1 - u), 
    
    which is on the negative x axis at x = mu - 1
    
    The final potential, 
    
    f = -G M_1 / r_1 - G M_2 / r_2 - omega^2 (x^2 + y^2)/2
    
    is written in these reduced coordinates and is defined as f =
    -V, defined in the function Vf() below.  This is written such
    that V is always positive.
    
    the positions of the Lagrange points are found iteratively using
    dV/dx and d^2 V/dx^2, also defined as functions below.
    
    therefore on the plots produced here, M_2 (the less massive star)
    is on the left always.
    """


    def __init__(self, mu, N):
        """ 
        Define an equipotential problem.  
        
        mu is the mass parameter

        N is the number of points for our grid

        """

        self.mu = mu
        self.N = N

        self.xmin = -2.0
        self.xmax = 2.0
        self.ymin = -2.0
        self.ymax = 2.0

        self.EPS = 1.e-12
        self.NITER = 100

        self.x = np.linspace(self.xmin, self.xmax, self.N, dtype=np.float64)
        self.y = np.linspace(self.ymin, self.ymax, self.N, dtype=np.float64)

        self.X, self.Y = np.meshgrid(self.x, self.y)

        # this is the dimensionless potential, in the z=0 plane
        self.V = self.Vf(self.X, self.Y)

    def get_L1(self):
        # find L1
        x0 = 0.0

        x_L1 = self._solve(-1.0, 1.0, x0)
        y_L1 = 0.0
        V_L1 = self.Vf(x_L1, y_L1)

        return x_L1, y_L1, V_L1

    def get_L2(self):
        # find L2
        x0 = -1.0

        x_L2 = self._solve(-1.0, -1.0, x0)
        y_L2 = 0.0
        V_L2 = self.Vf(x_L2, y_L2)

        return x_L2, y_L2, V_L2

    def get_L3(self):
        # find L3
        x0 = 1.0

        x_L3 = self._solve(1.0, 1.0, x0)
        y_L3 = 0.0
        V_L3 = self.Vf(x_L3, y_L3)

        return x_L3, y_L3, V_L3

    def get_L4(self):
        # L4
        x_L4 = self.mu - 0.5
        y_L4 = math.sin(math.pi/3.0)
        
        return x_L4, y_L4

    def get_L5(self):
        # L5
        x_L5 = self.mu - 0.5
        y_L5 = -math.sin(math.pi/3.0)

        return x_L5, y_L5

    def _solve(self, a, b, x0):

        for n in range(self.NITER):
            dVX = self.dVXdx(a, b, x0)
            d2VX = self.d2VXdx2(a, b, x0)

            x1 = x0 - dVX/d2VX
            err = abs(x1 - x0)/abs(x0 + self.EPS)

            if err < self.EPS: break
            x0 = x1

        return x0

    def Vf(self, x, y):
        V = (1.0 - self.mu)/np.sqrt( (x - self.mu)**2 + y**2 ) + \
            self.mu/np.sqrt( (x + 1.0 - self.mu)**2 + y**2 ) + \
            0.5*(x**2 + y**2)    
        return V

    def dVXdx(self, h1, h2, x):
        dVX = -h1*(1.0 - self.mu)/(x - self.mu)**2 - h2*self.mu/(x + 1.0 - self.mu)**2 + x
        return dVX

    def d2VXdx2(self, h1, h2, x):
        d2VX = 2.0*h1*(1.0 - self.mu)/(x - self.mu)**3 + 2.0*h2*self.mu/(x + 1.0 - self.mu)**3 + 1.0
        return d2VX


def make_plot(mu):

    # do imshow
    plt.clf()

    eq = Equipotentials(mu, 1024)

    print mu,  eq.V.min(), eq.V.max()
    plt.imshow(np.log10(eq.V), origin="lower", cmap="Accent",
               extent=[eq.xmin, eq.xmax, eq.ymin, eq.ymax])

    # draw contours -- these values seem reasonable for a range of mu's
    Vmin = 1.5 
    Vmax = 1000.0  # np.max(V)
    nC = 25
    
    C = np.logspace(math.log10(Vmin), math.log10(Vmax), nC)

    plt.contour(eq.x, eq.y, eq.V, C, colors="b")

    x_L1, y_L1, V_L1 = eq.get_L1()
    x_L2, y_L2, V_L2 = eq.get_L2()
    x_L3, y_L3, V_L3 = eq.get_L3()

    # special contours right through the lagrange points
    plt.contour(eq.x, eq.y, eq.V, [V_L1], colors="b")
    plt.contour(eq.x, eq.y, eq.V, [V_L2], colors="b")
    plt.contour(eq.x, eq.y, eq.V, [V_L3], colors="b")
 
    
    # mark the Lagrange points and write the names
    xeps = 0.025

    plt.scatter([x_L1], [y_L1], marker="x", color="r", s=50)
    plt.text(x_L1+xeps, y_L1+xeps, "L1", color="r")

    plt.scatter([x_L2], [y_L2], marker="x", color="r", s=50)
    plt.text(x_L2+xeps, y_L2+xeps, "L2", color="r")

    plt.scatter([x_L3], [y_L3], marker="x", color="r", s=50)
    plt.text(x_L3+xeps, y_L3+xeps, "L3", color="r")

    x_L4, y_L4 = eq.get_L4()
    plt.scatter([x_L4], [y_L4], marker="x", color="r", s=50)
    plt.text(x_L4+xeps, y_L4+xeps, "L4", color="r")

    x_L5, y_L5 = eq.get_L5()
    plt.scatter([x_L5], [y_L5], marker="x", color="r", s=50)
    plt.text(x_L5+xeps, y_L5+xeps, "L5", color="r")
       
    plt.axis([eq.xmin, eq.xmax, eq.ymin, eq.ymax])

    plt.title(r"Equipotentials, $\mu = M_2/(M_1 + M_2) = {:5.3f}$".format(mu) ,fontsize=12)
    plt.xlabel("$x/(a + b)$")
    plt.ylabel("$y/(a + b)$")

    f = plt.gcf()
    f.set_size_inches(10.8,10.8)

    plt.tight_layout()

    plt.savefig("equipotentials_mu_{:5.3f}.png".format(mu))




if __name__== "__main__":

    for mu in np.linspace(0.005, 0.5, 200):
        make_plot(mu)


