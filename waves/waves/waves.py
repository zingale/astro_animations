#!/bin/env python

import math
import numpy
import pylab

# illustrate the difference between wavelength and frequency

# M. Zingale 


def waves():

    tmax = 4.0
    dt = tmax/1200

    lambda_1 = 1.0
    lambda_2 = 0.25

    t = 0.0
    v = 2.0

    xmin = 0
    xmax = 3.0
    npts = 500

    x = numpy.arange(npts)*(xmax - xmin)/npts

    iframe = 0

    while (t < tmax):
         
        pylab.clf()

        pylab.subplot(211)

        y = numpy.sin(2*math.pi*(x - v*t)/lambda_1)

        xpt = 2.0
        ypt = math.sin(2*math.pi*(xpt - v*t)/lambda_1)
        
        pylab.plot(x,y,color="b")
        pylab.scatter([xpt],[ypt],color="b")

        # draw and annotate a dimension line
        pylab.annotate("", (xmin, 1.2), (xmin+lambda_1, 1.2),
                       arrowprops=dict(arrowstyle="|-|",mutation_scale=5.0))

        pylab.text(0.5*(xmin + xmin + lambda_1), 1.4, 
                   r"$\lambda = %3.2f\, \rm{cm}$" % (lambda_1), 
                   horizontalalignment='center')

        pylab.plot([2.0,2.0], [-1.5,1.5], "k--")

        pylab.axis([xmin,xmax,-1.6,1.6])
        pylab.axis("off")

        pylab.subplot(212)

        y = numpy.sin(2*math.pi*(x - v*t)/lambda_2)
        
        xpt = 2.0
        ypt = math.sin(2*math.pi*(xpt - v*t)/lambda_2)

        pylab.plot(x,y,color="r")
        pylab.scatter([xpt],[ypt],color="r")

        # draw and annotate a dimension line
        pylab.annotate("", (xmin, 1.2), (xmin+lambda_2, 1.2),
                       arrowprops=dict(arrowstyle="|-|",mutation_scale=5.0))

        pylab.text(0.5*(xmin + xmin + lambda_2), 1.4, 
                   r"$\lambda = %3.2f\, \rm{cm}$" % (lambda_2), 
                   horizontalalignment='center')

        pylab.plot([2.0,2.0], [-1.5,1.5], "k--")

        # draw the velocity vector 
        pylab.arrow(0.75*xmax, -1.5, 0.13*xmax, 0.0, color="k",
                    length_includes_head=True,
                    head_width = 0.15, head_length = 0.15, width=0.05)

        pylab.text( 0.90*xmax, -1.5, 
                    r"$v = %3.2f\, \rm{cm/s}$" % (v), 
                    verticalalignment='center')

        pylab.axis([xmin,xmax,-1.6,1.6])
        pylab.axis("off")

        pylab.text(0.01,-1.6, r"$t = %5.3f\, \rm{s}$" % (t))

        f = pylab.gcf()
        f.set_size_inches(7.2,7.2)

        outfile = "wave_%04d.png" % iframe
        pylab.savefig(outfile)

        iframe += 1
        t += dt



if __name__== "__main__":
    waves()


    
        
