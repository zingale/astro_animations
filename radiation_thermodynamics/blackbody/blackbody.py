import math
import numpy
import pylab
from matplotlib.ticker import MultipleLocator, NullLocator
# draw blackbody curves of varying temperature

# M. Zingale (2008-09-20)

# physical constants (cgs)
h = 6.626068e-27
c = 2.99792458e10
k = 1.3806503e-16

SMALL = 1.e-8


def blackbody():

    # define frequncy
    npts = 250

    fmin = 1.e4
    fmax = 1.e24

    dlogf = (math.log10(fmax) - math.log10(fmin))/(npts-1)

    f = 10.0**(numpy.arange(npts, dtype=numpy.float64)*dlogf + math.log10(fmin))


    npts_T = 500
    Tmin = 1.0
    Tmax = 1.e10
    dlogT = (math.log10(Tmax) - math.log10(Tmin))/(npts_T-1)

    T = 10.0**(numpy.arange(npts_T, dtype=numpy.float64)*dlogT + math.log10(Tmin))

    T_1 = 1.e2
    I_1 = (2.0*(h*f)*(f/c)*(f/c))/(numpy.exp(h*f/(k*T_1)) - 1.0)

    T_2 = 1.e4
    I_2 = (2.0*(h*f)*(f/c)*(f/c))/(numpy.exp(h*f/(k*T_2)) - 1.0)

    T_3 = 1.e6
    I_3 = (2.0*(h*f)*(f/c)*(f/c))/(numpy.exp(h*f/(k*T_3)) - 1.0)

    T_4 = 1.e8
    I_4 = (2.0*(h*f)*(f/c)*(f/c))/(numpy.exp(h*f/(k*T_4)) - 1.0)

    print T

    n = 0
    while (n < npts_T):
        print n, T[n]

        x = h*f/(k*T[n])
        I = (2.0*(h*f)*(f/c)*(f/c))/(numpy.exp(x) - 1.0)


        # plotting
        pylab.subplots_adjust(left=0.125,right=0.75,
                              bottom=0.1,top=0.9,wspace=0.1)       
        
        # turn on interactive mode 
        #pylab.ion()

        pylab.clf()
    
        pylab.loglog(f,I,color="g",linewidth=2.0)


        # reference lines
        if (T[n] >= T_1):
            pylab.loglog(f,I_1,color="#888888",linewidth=2.0, alpha=0.5)
            # use Wien's law to find the max
            f1max = c*T_1/0.29
            I1max = (2.0*h*(f1max/c)*(f1max/c)*f1max)/ \
                (numpy.exp(h*f1max/(k*T_1)) - 1.0)
            exp = int(math.log10(T_1))
            pre = T_1/10.0**exp
            if (pre == 1.0):
                pylab.text(1.1*f1max,1.1*I1max, 
                           r"$T = 10^{%d} \, \mathrm{K} $" % (exp), 
                           color="#666666")
            else:
                pylab.text(1.1*f1max,1.1*I1max, 
                           r"$T = %3.2f \times 10^{%d} \, \mathrm{K} $" % (pre, exp), 
                           color="#666666")

        if (T[n] >= T_2):
            pylab.loglog(f,I_2,color="#666666",linewidth=2.0, alpha=0.5)
            # use Wien's law to find the max
            f2max = c*T_2/0.29
            I2max = (2.0*h*(f2max/c)*(f2max/c)*f2max)/ \
                (numpy.exp(h*f2max/(k*T_2)) - 1.0)
            exp = int(math.log10(T_2))
            pre = T_2/10.0**exp
            if (pre == 1.0):
                pylab.text(1.1*f2max,1.1*I2max, 
                           r'$T = 10^{%d} \, \mathrm{K} $' % (exp), 
                           color="#666666")
            else:
                pylab.text(1.1*f2max,1.1*I2max, 
                           r"$T = %3.2f \times 10^{%d} \, \mathrm{K} $" % (pre, exp), 
                           color="#666666")


        if (T[n] >= T_3):
            pylab.loglog(f,I_3,color="#666666",linewidth=2.0, alpha=0.5)
            # use Wien's law to find the max
            f3max = c*T_3/0.29
            I3max = (2.0*h*(f3max/c)*(f3max/c)*f3max)/ \
                (numpy.exp(h*f3max/(k*T_3)) - 1.0)
            exp = int(math.log10(T_3))
            pre = T_3/10.0**exp
            if (pre == 1.0):
                pylab.text(1.1*f3max,1.1*I3max, 
                           r"$T = 10^{%d} \, \mathrm{K} $" % (exp), 
                           color="#666666")
            else:
                pylab.text(1.1*f3max,1.1*I3max, 
                           r"$T = %3.2f \times 10^{%d} \, \mathrm{K} $" % (pre, exp), 
                           color="#666666")


        if (T[n] >= T_4):
            pylab.loglog(f,I_4,color="#666666",linewidth=2.0, alpha=0.5)
            # use Wien's law to find the max
            f4max = c*T_4/0.29
            I4max = (2.0*h*(f4max/c)*(f4max/c)*f4max)/ \
                (numpy.exp(h*f4max/(k*T_4)) - 1.0)
            exp = int(math.log10(T_4))
            pre = T_4/10.0**exp
            if (pre == 1.0):
                pylab.text(1.1*f4max,1.1*I4max, 
                           r"$T = 10^{%d} \, \mathrm{K} $" % (exp), 
                           color="#666666")
            else:
                pylab.text(1.1*f4max,1.1*I4max, 
                           r"$T = %3.2f \times 10^{%d} \, \mathrm{K} $" % (pre, exp), 
                           color="#666666")




        # visible frequencies: 7.5e14 Hz (blue) to 4.2857e14 Hz (red)
        Imin = 1.e-18
        Imax = 1.e15

        pylab.fill([4.2857e14,4.2857e14,7.5e14,7.5e14,4.2857e14],
                   [Imin,Imax,Imax,Imin,Imin],alpha=0.20,facecolor="b")


        pylab.axis([fmin,fmax,Imin,Imax])
        pylab.xlabel("frequency [Hz]")
        pylab.ylabel("Intensity")
        pylab.title("Blackbody Radiation")



        # draw a thermometer 
        pylab.axes([0.885,0.1,0.015,0.8])

        # we do a lot of thin lines here, instead of using a thick
        # line width, because a thick line extends vertically as
        # well as horizontally, and does not line up to the right place
        # on the axis
        therm_xmax = 0.01
        iter = 0
        while (iter < 10):
            pylab.semilogy([iter*therm_xmax/10.0,iter*therm_xmax/10.0],
                           [T[0],T[n]], color='r',linewidth=2)
            iter += 1


        pylab.axis([0,therm_xmax,T[0],T[npts_T-1]])

        axis = pylab.gca()
        #majorLocator = MultipleLocator(20)
        majorLocator = NullLocator()
        axis.xaxis.set_major_locator(majorLocator)

        pylab.ylabel("temperature [K]")

        fig = pylab.gcf()
        fig.set_size_inches(7.5,6.0)

        outfile = "blackbody_%03d.png" % n
        pylab.savefig(outfile)

        n += 1



if __name__== "__main__":
    blackbody()


    
        
