import math
import numpy
import pylab
import random

# simulate radioactive decay by populating a grid with a number of
# markers and having each marker have a 50/50 chance of "decaying"
# each half-life

class marker:

    def __init__(self, xc, yc):

        # a marker is indicated by its center (xc,yc).
        self.xc = xc
        self.yc = yc
        
        # the state indicates whether it is the original object (1)
        # or decayed into its daughter product (0)
        self.state = 1
        

    def decay(self):
    
        if (self.state == 1):
            
            # random.random() returns a number in the range [0.0, 1.0)
            if (random.random() >= 0.5):
                self.state = 0


def num_decayed(markers):

    nd = 0
    n = 0
    while (n < len(markers)):
        if (markers[n].state == 0):
            nd += 1

        n += 1

    return nd


def radioactive_decay():

    # define the number of markers in x and y
    nx = 50
    ny = 50

    # estimate the number of half-lifes needed to decay all markers
    nest = int(math.log(nx*ny)/math.log(2))

    # allocate storage for the number of markers that are decayed at
    # each half-life.  For safety, we allocate space for 2x our estimate
    hist_decay = numpy.zeros(2*nest)
    

    # define the length of a marker side
    L = 0.8

    # create a list of marker objects, one at each grid location
    markers = []
    i = 0
    while (i < nx):

        j = 0
        while (j < ny):

            markers.append(marker(i, j))

            j += 1

        i += 1


    # loop over half-lives, and re-evaluate the marker state
    t = 0
    while (t < 2*nest):

        pylab.clf()

        # the margins are funny -- we pick them to ensure that the
        # plot size is an integer multiple of the number of markers in
        # each dimension
        pylab.subplots_adjust(left=0.0493333,  right=0.9506666,
                              bottom=0.0493333,top=0.9506666)

        # draw the current state
        n = 0
        while (n < len(markers)):

            if (markers[n].state == 1):
                c = "r"
            else:
                c = "w"

            pylab.fill([markers[n].xc-L/2, markers[n].xc-L/2,
                        markers[n].xc+L/2, markers[n].xc+L/2,
                        markers[n].xc-L/2],
                       [markers[n].yc-L/2, markers[n].yc+L/2,
                        markers[n].yc+L/2, markers[n].yc-L/2,
                        markers[n].yc-L/2], 
                       c)

            n += 1

            
        nd = num_decayed(markers)
        hist_decay[t] = nd

        ax = pylab.axis([-1,nx+1,-1,ny+1])
        pylab.axis("off")

        pylab.text(-1,-1, "number decayed = %d" % (nd),horizontalalignment="left", verticalalignment="top")

        f = pylab.gcf()
        f.set_size_inches(7.5,7.5)

        outfile = "radioactive_decay_%04d.png" % t
        pylab.savefig(outfile)

        print t, nd

        # if all are decayed, stop making plots
        if (nd == nx*ny):
            break

        # now give each marker a chance to "decay"
        n = 0
        while (n < len(markers)):
            markers[n].decay()
            n += 1
        
        t += 1


    # now make a plot of the number that decayed as a function of half-life
    pylab.clf()
    pylab.subplots_adjust(left=0.1,   right=0.9,
                          bottom=0.1, top=0.9)


    tsmooth = numpy.arange(1000)*t/1000.
    pred = nx*ny*(1.0/2.0)**tsmooth

    pylab.plot(tsmooth, pred, color="b")

    time = numpy.arange(t-1)
    pylab.scatter(time, nx*ny-hist_decay[0:t-1], color="b", label="parent")
    pylab.scatter(time, hist_decay[0:t-1], color="r", label="daughter")

    pylab.axis([0,numpy.max(time)+1,0,1.1*nx*ny])
    pylab.xlabel("half life")

    leg = pylab.legend(loc=2)
    ltext = leg.get_texts()
    pylab.setp(ltext, fontsize='small')
    leg.draw_frame(0)



    outfile = "radioactive_decay_%04d.png" % (t+1)
    pylab.savefig(outfile)



if __name__== "__main__":
    radioactive_decay()


    
        
