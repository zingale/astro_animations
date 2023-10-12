import numpy as np
import pylab

# illustrate the wind-up problem of spiral arms

# M. Zingale (2009-04-20)


# define the rotation curve 
def vel(r):
    ar = np.abs(r)
    return np.piecewise(ar, [ar<=5.0, ar>5.0 and ar<6.0],
                        [lambda x: 1.4*x*np.exp(-x/20.0),
                         lambda x: (5.0 - 7.0*np.exp(-0.25))*(x - 5.0) + 7.0*np.exp(-0.25),
                         lambda x: 5])


def spiral():

    rmin = -20
    rmax = 20
    npts = 1000

    r = np.arange(npts, dtype=np.float64)*(rmax - rmin)/(npts - 1.0) + rmin

    omega = np.zeros( npts )

    n = 0
    while (n < npts):
        omega[n] = np.abs(vel(r[n])/r[n])
        n += 1

    t = np.arange(0.0, 20.0, 0.04)
    t = 0.0
    tmax = 20.0
    dt = 0.04

    iframe = 0

    while (t < tmax):
 
        
        pylab.clf()

        x = r*np.cos(omega*t)
        y = r*np.sin(omega*t)

        pylab.plot(x,y,color="b")


        pylab.axis([1.2*rmin,1.2*rmax,1.2*rmin,1.2*rmax])
        pylab.axis("off")

        f = pylab.gcf()
        f.set_size_inches(6.0,6.0)

        outfile = "spiral_%04d.png" % iframe
        pylab.savefig(outfile)

        iframe += 1
        t += dt



if __name__== "__main__":
    spiral()
