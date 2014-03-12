import math
import numpy
import pylab


# the wavefront class will hold the location and time that a wavefront
# was emitted
class wavefront:

    def __init__ (self, x_emit, y_emit, w, t_emit):
        self.x_emit = x_emit
        self.y_emit = y_emit
        self.w      = w           # wave propagation speed
        self.t_emit = t_emit



def doppler():

    # emitter velocity (in x-direction)
    vel1 = 1.0
    vel2 = 0.5

    # emitter initial coords
    x_init = 0.0
    y_init = 0.0

    # wave velocity
    w = 2.0

    # wave frequency (# of peaks per second)
    f = 3.0

    # maximum time
    tmax = 10.0
    dt = 0.01

    # create a list of wavefront objects that we can refer to when we
    # want to plot things.  There are f wavefronts emitted per second,
    # so the total number of wavefronts is tmax*f
    t = 0

    wavefronts1 = []
    while (t <= tmax):
        
        x_emit = x_init + vel1*t
        y_emit = y_init

        wavefronts1.append(wavefront(x_emit, y_emit, w, t))
        
        t += 1/f


    t = 0

    wavefronts2 = []
    while (t <= tmax):
        
        x_emit = x_init + vel2*t
        y_emit = y_init

        wavefronts2.append(wavefront(x_emit, y_emit, w, t))
        
        t += 1/f

        


    xmax = x_init + max(vel1,vel2)*tmax

    # we will be drawing circles, so make an array with the polar angle
    npts = 360
    theta = numpy.arange(npts)*2*math.pi/(npts-1)


    # step forward in time (by dt) and draw any wavefronts that have
    # been emitted
    iframe = 0
    t = 0
    while (t <= tmax):

        pylab.clf()

        pylab.subplot(211)

        x_source = x_init + vel1*t
        y_source = y_init

        # plot the sources's path
        pylab.plot([-1.2*xmax,1.2*xmax],[y_init,y_init],'k:')

        # draw the source
        pylab.scatter([x_source],[y_source], color='b')

        # loop over the wavefronts, and draw any that have been
        # emitted so far
        n = 0
        while (n < len(wavefronts1)):

            if (wavefronts1[n].t_emit > t):
                break

            r_front = wavefronts1[n].w*(t - wavefronts1[n].t_emit)

            # wavefronts are circles centered on their emitted coordinates
            x_front = wavefronts1[n].x_emit + r_front*numpy.cos(theta)
            y_front = wavefronts1[n].y_emit + r_front*numpy.sin(theta)
            
            pylab.plot(x_front, y_front, color='r')
            
            n += 1


        pylab.plot([-1.2*xmax,1.2*xmax],[-0.8*xmax,-0.8*xmax], color='k',lw=2)
        pylab.axis([-1.2*xmax,1.2*xmax,-0.8*xmax,0.8*xmax])
        pylab.axis("off")


        pylab.subplot(212)

        x_source = x_init + vel2*t
        y_source = y_init

        # plot the sources's path
        pylab.plot([-1.2*xmax,1.2*xmax],[y_init,y_init],'k:')

        # draw the source
        pylab.scatter([x_source],[y_source], color='b')

        # loop over the wavefronts, and draw any that have been
        # emitted so far
        n = 0
        while (n < len(wavefronts2)):

            if (wavefronts2[n].t_emit > t):
                break

            r_front = wavefronts2[n].w*(t - wavefronts2[n].t_emit)

            # wavefronts are circles centered on their emitted coordinates
            x_front = wavefronts2[n].x_emit + r_front*numpy.cos(theta)
            y_front = wavefronts2[n].y_emit + r_front*numpy.sin(theta)
            
            pylab.plot(x_front, y_front, color='g')
            
            n += 1
            

        pylab.plot([-1.2*xmax,1.2*xmax],[0.8*xmax,0.8*xmax], color='k',lw=2)
        pylab.axis([-1.2*xmax,1.2*xmax,-0.8*xmax,0.8*xmax])
        pylab.axis("off")

        
        pylab.subplots_adjust(left=0,right=1.0,bottom=0,top=1.0,
                              wspace=0,hspace=0)


        f = pylab.gcf()
        f.set_size_inches(4.5,6.0)

        outfile = "doppler_%04d.png" % iframe
        pylab.savefig(outfile)



        t += dt
        iframe += 1
    

    


if __name__== "__main__":
    doppler()

