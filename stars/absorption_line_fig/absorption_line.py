import math
import numpy
import pylab

def absorption_line():

    ax = pylab.subplot(111)    

    npts = 200
    lambda_x = numpy.arange(npts, dtype=numpy.float64)/npts
    lambda_0 = 0.5
    W = 0.075
    A = 0.9

    flux = 1.0 - A*numpy.exp(-(lambda_x - lambda_0)**2/(2*W**2))

    # if the gaussian is a exp{-(x-b)^2 / 2c^2 }, the area is ac sqrt{2 pi}
    area = A*W*math.sqrt(2.0*math.pi)

    width = area/A

    print lambda_x
    print flux

    pylab.plot(lambda_x, flux, 'r')

    pylab.plot([lambda_0 - width/2, lambda_0 - width/2], [0,1], 'b--')
    pylab.plot([lambda_0 + width/2, lambda_0 + width/2], [0,1], 'b--')

    pylab.ylabel("flux")

    tlabels = [r"$\lambda_0$"]
    tpos = [lambda_0]

    locs, labels = pylab.xticks(tpos, tlabels)

    f = pylab.gcf()
    f.set_size_inches(6.0,6.0)

    pylab.savefig("absorption_line.png")



if __name__== "__main__":
    absorption_line()

