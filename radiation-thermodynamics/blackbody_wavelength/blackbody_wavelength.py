import math
import numpy
import pylab
from matplotlib.ticker import MultipleLocator, NullLocator

# draw blackbody curves of varying temperature
# here we plot I_lambda

# M. Zingale (2008-09-20)

# physical constants (cgs)
h = 6.63e-27
c = 3.e10
k = 1.38e-16



def blackbody():

    # define wavelength
    npts = 250

    lmin = 1e-14
    lmax = 1e5

    dlogl = (math.log10(lmax) - math.log10(lmin))/(npts-1)

    l = 10.0**(numpy.arange(npts)*dlogl + math.log10(lmin))


    npts_T = 500
    Tmin = 1.0
    Tmax = 1.e10
    dlogT = (math.log10(Tmax) - math.log10(Tmin))/(npts_T-1)

    T = 10.0**(numpy.arange(npts_T)*dlogT + math.log10(Tmin))

    T_1 = 1.e2
    I_1 = (2.0*h*c*c/(l**5))/(numpy.exp(h*c/(l*k*T_1)) - 1.0)

    T_2 = 1.e4
    I_2 = (2.0*h*c*c/(l**5))/(numpy.exp(h*c/(l*k*T_2)) - 1.0)

    T_3 = 1.e6
    I_3 = (2.0*h*c*c/(l**5))/(numpy.exp(h*c/(l*k*T_3)) - 1.0)

    T_4 = 1.e8
    I_4 = (2.0*h*c*c/(l**5))/(numpy.exp(h*c/(l*k*T_4)) - 1.0)

    print T

    n = 0
    while (n < npts_T):
        print n, T[n]

        I = (2.0*h*c*c/(l**5))/(numpy.exp(h*c/(l*k*T[n])) - 1.0)


        # plotting
        pylab.subplots_adjust(left=0.125,right=0.75,
                              bottom=0.1,top=0.9,wspace=0.1)       
        
        # turn on interactive mode 
        #pylab.ion()

        pylab.clf()
    
        pylab.loglog(l,I,color="g",linewidth=2.0)


        # reference lines
        if (T[n] >= T_1):
            pylab.loglog(l,I_1,color="#888888",linewidth=2.0, alpha=0.5)
            # use Wien's law to find the max
            l1max = 0.29/T_1
            I1max = (2.0*h*c*c/(l1max**5))/(numpy.exp(h*c/(l1max*k*T_1)) - 1.0)
            exp = int(math.log10(T_1))
            pre = T_1/10.0**exp
            if (pre == 1.0):
                pylab.text(1.1*l1max,1.1*I1max, 
                           r"$T = 10^{%d} \, \mathrm{K} $" % (exp), 
                           color="#666666")
            else:
                pylab.text(1.1*l1max,1.1*I1max, 
                           r"$T = %3.2f \times 10^{%d} \, \mathrm{K} $" % (pre, exp), 
                           color="#666666")

        if (T[n] >= T_2):
            pylab.loglog(l,I_2,color="#888888",linewidth=2.0, alpha=0.5)
            # use Wien's law to find the max
            l2max = 0.29/T_2
            I2max = (2.0*h*c*c/(l2max**5))/(numpy.exp(h*c/(l2max*k*T_2)) - 1.0)
            exp = int(math.log10(T_2))
            pre = T_2/10.0**exp
            if (pre == 1.0):
                pylab.text(1.1*l2max,1.1*I2max, 
                           r"$T = 10^{%d} \, \mathrm{K} $" % (exp), 
                           color="#666666")
            else:
                pylab.text(1.1*l2max,1.1*I2max, 
                           r"$T = %3.2f \times 10^{%d} \, \mathrm{K} $" % (pre, exp), 
                           color="#666666")


        if (T[n] >= T_3):
            pylab.loglog(l,I_3,color="#888888",linewidth=2.0, alpha=0.5)
            # use Wien's law to find the max
            l3max = 0.29/T_3
            I3max = (2.0*h*c*c/(l3max**5))/(numpy.exp(h*c/(l3max*k*T_3)) - 1.0)
            exp = int(math.log10(T_3))
            pre = T_3/10.0**exp
            if (pre == 1.0):
                pylab.text(1.1*l3max,1.1*I3max, 
                           r"$T = 10^{%d} \, \mathrm{K} $" % (exp), 
                           color="#666666")
            else:
                pylab.text(1.1*l3max,1.1*I3max, 
                           r"$T = %3.2f \times 10^{%d} \, \mathrm{K} $" % (pre, exp), 
                           color="#666666")


        if (T[n] >= T_4):
            pylab.loglog(l,I_4,color="#888888",linewidth=2.0, alpha=0.5)
            # use Wien's law to find the max
            l4max = 0.29/T_4
            I4max = (2.0*h*c*c/(l4max**5))/(numpy.exp(h*c/(l4max*k*T_4)) - 1.0)
            exp = int(math.log10(T_4))
            pre = T_4/10.0**exp
            if (pre == 1.0):
                pylab.text(1.1*l4max,1.1*I4max, 
                           r"$T = 10^{%d} \, \mathrm{K} $" % (exp), 
                           color="#666666")
            else:
                pylab.text(1.1*l4max,1.1*I4max, 
                           r"$T = %3.2f \times 10^{%d} \, \mathrm{K} $" % (pre, exp), 
                           color="#666666")




        # visible wavelengths: 400 nm to 700 nm  (4e-5 cm to 7e-5 cm)
        Imin = 1.e-6
        Imax = 1.e48

        pylab.fill([4.e-5,4.e-5,7.e-5,7.e-5,4.e-5],
                   [Imin,Imax,Imax,Imin,Imin],alpha=0.20,facecolor="b")


        pylab.axis([lmin,lmax,Imin,Imax])
        pylab.xlabel("wavelength [cm]")
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


    
        
