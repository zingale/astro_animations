import math
import numpy
import pylab
import random

def random_walk():

    # length of single step
    l = 1.0

    # maximum radius of domain
    R = 25.0

    # maximum number of steps
    N = 5000

    # seed for the random number generator
    seed = 1000
    random.seed(seed)

    # take steps, draw a segment in a random direction, and save the frame
    pylab.clf()

    # draw a circle to indicate the extent of the domain
    npts = 360
    theta = numpy.arange(npts)*2*math.pi/(npts-1)

    pylab.plot(R*numpy.cos(theta), R*numpy.sin(theta), color='k')

    pylab.subplots_adjust(left=0,right=1.0,bottom=0,top=1.0)

    # set the inital position to be the origin
    x_0 = 0.0
    y_0 = 0.0

    n = 0
    while (n < N):

        # compute a random angle
        angle = 2.0*math.pi*random.random()

        # compute the end coordinates of the segment
        x_1 = l*math.cos(angle) + x_0
        y_1 = l*math.sin(angle) + y_0
                
        pylab.plot([x_0,x_1], [y_0,y_1], color='r')

        pylab.axis([-1.1*R,1.1*R,-1.1*R,1.1*R])
        pylab.axis("off")

        f = pylab.gcf()
        f.set_size_inches(5.0,5.0)

        outfile = "random_walk_%04d.png" % n
        pylab.savefig(outfile)

        # have we hit the edge of our domain?
        if (math.sqrt(x_1**2 + y_1**2) >= R):
            break
        else:
            x_0 = x_1
            y_0 = y_1

        n += 1

        

if __name__== "__main__":
    random_walk()

