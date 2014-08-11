import math
import numpy
import pylab

# illustrate the difference between wavelength and frequency

# M. Zingale (2009-01-28)


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
        pylab.plot([xmin,xmin],[1.1,1.3],color='k')
        pylab.plot([xmin+lambda_1,xmin+lambda_1],[1.1,1.3],color='k')
        pylab.plot([xmin,xmin+lambda_1],[1.2,1.2],color='k')
        pylab.text(0.5*(xmin + xmin + lambda_1), 1.4, r"$\lambda = %3.2f\, \rm{cm}$" % (lambda_1), horizontalalignment='center')

        pylab.plot([2.0,2.0],[-1.5,1.5],"k--")

        pylab.axis([xmin,xmax,-1.6,1.6])
        pylab.axis("off")


        pylab.subplot(212)

        y = numpy.sin(2*math.pi*(x - v*t)/lambda_2)
        
        xpt = 2.0
        ypt = math.sin(2*math.pi*(xpt - v*t)/lambda_2)

        pylab.plot(x,y,color="r")
        pylab.scatter([xpt],[ypt],color="r")

        # draw and annotate a dimension line
        pylab.plot([xmin,xmin],[1.1,1.3],color='k')
        pylab.plot([xmin+lambda_2,xmin+lambda_2],[1.1,1.3],color='k')
        pylab.plot([xmin,xmin+lambda_2],[1.2,1.2],color='k')
        pylab.text(0.5*(xmin + xmin + lambda_2), 1.4, r"$\lambda = %3.2f\, \rm{cm}$" % (lambda_2), horizontalalignment='center')

        pylab.plot([2.0,2.0],[-1.5,1.5],"k--")

        # draw the velocity vector -- I should be able to do this with 
        # the pylab annotate function, but I couldn't get it to work.
        pylab.plot( [0.75*xmax,0.88*xmax], [-1.5,-1.5], color='k')
        pylab.plot( [0.845*xmax,0.88*xmax], [-1.6,-1.5], color='k')
        pylab.plot( [0.845*xmax,0.88*xmax], [-1.4,-1.5], color='k')
        pylab.text( 0.90*xmax, -1.5, r"$v = %3.2f\, \rm{cm/s}$" % (v), verticalalignment='center')

        pylab.axis([xmin,xmax,-1.6,1.6])
        pylab.axis("off")

        pylab.text(0.01,-1.6, r"$t = %5.3f\, \rm{s}$" % (t))

        f = pylab.gcf()
        f.set_size_inches(6.0,6.0)

        outfile = "wave_%04d.png" % iframe
        pylab.savefig(outfile)

        iframe += 1
        t += dt



if __name__== "__main__":
    waves()


    
        
