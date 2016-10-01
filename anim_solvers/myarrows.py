import anim_solvers.stick_figure as sf
import numpy as np
import matplotlib.pyplot as plt


class Cap(object):
    """ caps are simple two lines making a 'V' that we rotate """

    def __init__(self, c, angle, L=1.0):

        self.L = L
        self.c = c
        self.angle = angle

        self.cap_x = [-0.5*self.L, 0, 0.5*self.L]
        self.cap_y = [-self.L, 0, -self.L]

    def draw(self, color="k", ls="-"):

        # arrow at start
        xa = []
        ya = []
        for x, y in zip(self.cap_x, self.cap_y):
            p = sf._rotate((x, y), (0,0), self.angle)
            xa += [self.c[0] + p[0]]
            ya += [self.c[1] + p[1]]

        plt.plot(xa, ya, color=color, ls=ls)


class ArcArrow(object):
    """ draw an arc with arrows on the end, in data coordinates """

    def __init__(self, center, radius, theta_start=90, theta_end=270):
        self.x0 = center[0]
        self.y0 = center[1]

        self.R = radius

        self.theta0 = np.radians(theta_start)
        self.theta1 = np.radians(theta_end)

    def draw(self, color="k", ls="-"):

        theta = np.linspace(self.theta0, self.theta1, 100)

        plt.plot(self.x0 + self.R*np.cos(theta),
                 self.y0 + self.R*np.sin(theta), color=color, ls=ls)


        start_point = (self.x0 + self.R*np.cos(theta[0]),
                       self.y0 + self.R*np.sin(theta[0]))

        end_point = (self.x0 + self.R*np.cos(theta[-1]),
                     self.y0 + self.R*np.sin(theta[-1]))


        # arrow at start
        c = Cap(start_point, self.theta0+np.pi, L=0.05*self.R)
        c.draw(color=color, ls=ls)

        # arrow at end
        c = Cap(end_point, self.theta1, L=0.05*self.R)
        c.draw(color=color, ls=ls)


if __name__ == "__main__":
    a = ArcArrow((0,0), 1.0, theta_start=135, theta_end=270)
    a.draw()

    plt.axis("off")

    ax = plt.gca()
    ax.set_aspect("equal", "datalim")

    plt.axis([-1.2*a.R, 1.2*a.R, -1.2*a.R, 1.2*a.R])

    f = plt.gcf()
    f.set_size_inches(6.0, 6.0)

    plt.savefig("atest.png")


        
