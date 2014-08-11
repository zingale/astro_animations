import math
import numpy
import pylab

# show how ellipses are drawn

# M. Zingale (2008-09-02)
#
#      2009-02-14: fixed the aspect ratio

def ellipse():

    # theta ranges from 0 to 2pi
    npts = 250
    theta = numpy.arange(npts)*2.0*math.pi/(npts-1)

    e = 0.5
    a = 1.0

    n = 0
 
    r = a*(1.0 - e*e)/(1.0 + e*numpy.cos(theta))

    x = r*numpy.cos(theta)
    y = r*numpy.sin(theta)


    # plotting
    n = 0
    while (n < npts):

        # turn on interactive mode 
        #pylab.ion()

        pylab.clf()

        ax = pylab.gca()
        ax.set_aspect("equal", "datalim")
        
        # draw the ellipse
        pylab.plot(x,y,color="b")

        # draw our current point
        pylab.scatter([x[n]],[y[n]], color="k", s=75)

        # second foci
        pylab.scatter([-2.0*a*e], [0], color="g", marker="x", s=200)

        # primary foci
        pylab.scatter([0],           [0], color="r", marker="x", s=200)

        # draw lines connecting the foci to the current point
        pylab.plot([0,x[n]],[0,y[n]], color="r")
        pylab.plot([-2.0*a*e,x[n]],[0,y[n]], color="g")

        pylab.axis([-2.5,1.5,-2.,2.])

        len1 = math.sqrt( (x[n] - 0)**2 + (y[n] - 0)**2)
        len2 = math.sqrt( (x[n] - (-2.0*a*e))**2 + (y[n] - 0)**2)

        pylab.title("Ellipse, eccenticity = %5.3f" % e)

        pylab.text(-1.5,-1.25, "r length: %5.3f" % len1,color="r")
        pylab.text(-1.5,-1.5, "r' length: %5.3f" % len2, color="g")
        pylab.text(-1.5,-1.75, "r + r' = %5.3f" % (len1 + len2) )


        print n, len1, len2, len1+len2

        f = pylab.gcf()
        f.set_size_inches(6.0,6.0)

        outfile = "ellipsedraw_%03d.png" % n
        pylab.savefig(outfile)

        n += 1


if __name__== "__main__":
    ellipse()


    
        
