import math
import numpy
import pylab

# draw ellipses of varying eccentricity

# M. Zingale (2008-09-02)
#
#      2009-02-14: fix aspect ratio

def ellipse():

    # theta ranges from 0 to 2pi
    npts = 250
    theta = numpy.arange(npts)*2.0*math.pi/(npts-1)

    n_ecc = 100
    e_tmp = numpy.arange(n_ecc)*.95/n_ecc
    e = numpy.zeros(2*n_ecc)
    e[0:n_ecc] = e_tmp[:]
    e[n_ecc:] = e_tmp[::-1]

    a = 1.0

    n = 0
    while (n < 2*n_ecc):
 
        r = a*(1.0 - e[n]*e[n])/(1.0 + e[n]*numpy.cos(theta))

        x = r*numpy.cos(theta)
        y = r*numpy.sin(theta)


        # plotting
        
        # turn on interactive mode 
        #pylab.ion()

        pylab.clf()

        ax = pylab.gca()
        ax.set_aspect("equal", "datalim")

        pylab.plot(x,y,color="b")

        # second foci
        pylab.scatter([-2.0*a*e[n]], [0], color="r", marker="x", s=200)
        # primary foci
        pylab.scatter([0],[0],color="g",marker="x",s=200)



        pylab.axis([-2.5,1.5,-2.,2.])

        pylab.text(-1.5,-1.5,"a = %5.3f, e = %6.4f" % (a, e[n]))

        f = pylab.gcf()
        f.set_size_inches(6.0,6.0)

        outfile = "ellipse_%03d.png" % n
        pylab.savefig(outfile)

        n += 1



if __name__== "__main__":
    ellipse()


    
        
