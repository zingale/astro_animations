import math
import numpy
import pylab

# show the ellipse geometry

# M. Zingale (2009-02-13)


def ellipse():

    # theta ranges from 0 to 2pi
    npts = 250
    theta = numpy.arange(npts)*2.0*math.pi/(npts-1)

    e = 0.5
    a = 1.0
    b = a*math.sqrt(1-e*e)

    n = npts/5
 
    r = a*(1.0 - e*e)/(1.0 + e*numpy.cos(theta))

    x = r*numpy.cos(theta)
    y = r*numpy.sin(theta)


    # plotting
    pylab.clf()

    ax = pylab.gca()
    ax.set_aspect("equal", "datalim")
    pylab.axis("off")
        
    # Draw the ellipse
    pylab.plot(x,y,color="b")

    # draw our current point
    pylab.scatter([x[n]],[y[n]], color="k", s=25)

    # second foci
    pylab.scatter([-2.0*a*e],[0], color="g", marker="x", s=200)

    # primary foci
    pylab.scatter([0],       [0], color="g", marker="x", s=200)

    # center
    pylab.scatter([-a*e],       [0], color='k', marker="x", s=200)

    # draw the semi-major axis
    pylab.plot([-a*e,a*(1-e)],[0,0], color='r')
    pylab.text(0.5*a - a*e, -0.15*b, "a", color='r')

    # draw the semi-minor axis
    pylab.plot([-a*e,-a*e], [0,-b], color='k')
    pylab.text(-a*e-0.15*a*e, -0.5*b, "b", color='k')

    pylab.plot([-a*e, -2*a*e], [0,0], color="0.5")
    pylab.text(-1.5*a*e, -0.15*b, "ae", color="0.5")

    # draw lines connecting the foci to the current point
    pylab.plot([0,x[n]],[0,y[n]], color="g")
    pylab.text(0, 0.3*b, "r", color="g")

    pylab.plot([-2.0*a*e,x[n]],[0,y[n]], color="g")
    pylab.text(-a*e-0.1*a, 0.3*b, r"r$^\prime$", color="g")
    # indicate the angle
    pylab.text(0.075*a, 0.05*b, r"$\theta$", color="k")
    
    pylab.axis([-1.75,.75,-1.25,1.25])
 
    f = pylab.gcf()
    f.set_size_inches(6.0,6.0)

    outfile = "ellipse_geom.png"
    pylab.savefig(outfile)


if __name__== "__main__":
    ellipse()


    
        
