from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

def init_figure():
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.set_aspect("equal")
    ax.axis("off")
    return ax, fig


# rotation matrix:
# http://stackoverflow.com/questions/12148351/efficiently-rotate-a-set-of-points-with-a-rotation-matrix-in-numpy

# example:
# http://stackoverflow.com/questions/14824893/how-to-draw-diagrams-like-this


class Sphere:
    def __init__(self, R=1.0, center=(0., 0., 0.),
                 color="0.5", alpha=0.25, ax=None):

        if ax is None:
            ax, fig = init_figure()

        self.R = R
        self.center = center
        self.color = color
        self.alpha = alpha

        # see http://stackoverflow.com/questions/11140163/python-matplotlib-plotting-a-3d-cube-a-sphere-and-a-vector
        theta, phi = np.mgrid[0:np.pi:25j, 0:2.0*np.pi:50j]

        x = R*np.sin(theta)*np.cos(phi)
        y = R*np.sin(theta)*np.sin(phi)
        z = R*np.cos(theta)

        ax.plot_surface(x, y, z, color=color, alpha=alpha,
                        linewidth=0,
                        rstride=1, cstride=1)


class CPlane:
    """ a circular plane centered at c, with radius R, and normal n """
    def __init__(self, c=(0.,0.,0.), R=1.0, n=(0.,0.,1.),
                 color="0.5", alpha=0.75, ax=None):

        if ax == None:
            ax, fig = init_figure()


        # plane equation through x0, y0, z0 and normal n = (a, b, c)
        # a(x-x0) + b(y-y0) + c(z-z0) = 0

        if not n[2] == 0:
            r, phi = np.mgrid[0:R:20j, 0:2.0*np.pi:50j]

            x = r*np.cos(phi)
            y = r*np.sin(phi)

            z = c[2] - n[0]/n[2]*(x - c[0]) - n[1]/n[2]*(y - c[1])

        ax.plot_surface(x, y, z, color=color, alpha=alpha,
                        linewidth=0,
                        rstride=1, cstride=1, antialiased=False)


ax, fig = init_figure()

s = Sphere(ax=ax, color="b", alpha=0.1)
e = Sphere(ax=ax, color="b", alpha=1.0, R=0.1)

p = CPlane(ax=ax, n=(1, 0, 1))


plt.savefig("sphere.png")
