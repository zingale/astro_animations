import math
import numpy
import pylab

# illustrate the wind-up problem of spiral arms

# M. Zingale (2009-04-20)


# define the rotation curve 
def vel(r):

    ar = math.fabs(r)
    if (ar <= 5.0):
        v = 1.4*ar*math.exp(-ar/20.0)
    elif (ar > 5.0 and ar < 6.0):
        max = 1.4*5*math.exp(-5.0/20.0)
        min = 5.0
        m = (max - min)/(-1.0)
        v = m*(ar - 5.0) + max
    else:
        v = 5.0

    return v



def spiral():

    rmin = -20
    rmax = 20
    npts = 1000

    r = numpy.arange(npts, dtype=numpy.float64)*(rmax - rmin)/(npts - 1.0) + rmin

    omega = numpy.zeros( (npts) )

    n = 0
    while (n < npts):
        omega[n] = math.fabs(vel(r[n])/r[n])
        n += 1

        
    t = 0.0
    tmax = 20.0
    dt = 0.04

    iframe = 0

    while (t < tmax):
 
        
        pylab.clf()

        x = r*numpy.cos(omega*t)
        y = r*numpy.sin(omega*t)

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


    
        
