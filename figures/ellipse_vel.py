import math
import numpy
import pylab

# show the velocity vectors at perihelion and aphelion for ellipse geometry

# M. Zingale (2010-09-08)


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

    # primary foci
    pylab.scatter([0],       [0], color="k", marker="x", s=200)

    # draw arrow to perihelion
    eps = 0.033*a
    pylab.plot([0,a*(1-e)],[0,0], color='r')
    pylab.plot([a*(1-e),a*(1-e)-eps],[0,eps], color='r')
    pylab.plot([a*(1-e),a*(1-e)-eps],[0,-eps], color='r')
    pylab.text(0.5*a*(1-e), -0.15*b, "r", color='r')

    
    # draw velocity at perihelion
    vel = 0.66*a
    pylab.plot([a*(1-e),a*(1-e)], [0, vel], color='g')
    pylab.plot([a*(1-e),a*(1-e)-eps], [vel, vel-eps], color='g')
    pylab.plot([a*(1-e),a*(1-e)+eps], [vel, vel-eps], color='g')
    pylab.text(a*(1-e)+0.05*b,0.5*vel, "v", color='g')
    
    pylab.axis([-1.75,.75,-1.25,1.25])
 

    f = pylab.gcf()
    f.set_size_inches(6.0,6.0)

    outfile = "ellipse_vel.png"
    pylab.savefig(outfile)


if __name__== "__main__":
    ellipse()


    
        
