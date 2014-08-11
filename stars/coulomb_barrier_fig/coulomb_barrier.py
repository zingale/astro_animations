import math
import numpy
import pylab

def coulomb_barrier():

    ax = pylab.subplot(111)    
    
    npts = 200
    r = numpy.arange(npts, dtype=numpy.float64)/npts

    SMALL = 1.e-16

    # Coulomb potential
    C = 1/(r + SMALL)

    # interior, strong potential
    V = numpy.zeros(npts) - 2.0

    # merge the two
    U = numpy.where(C < 10.0, C, V)

    pylab.plot(r, U, 'r', linewidth=2)


    tlabels = [r"$\lambda_0$"]
    tpos = [0.1]

    # draw the axes manually
    pylab.axis("off")
    
    # y-axis
    pylab.arrow(0, 1.1*min(U), 0, 1.1*max(U) - 1.1*min(U), 
                fc='k', head_width=0.015, head_length=0.18)

    pylab.text(-0.05*max(r), 1.1*max(U), r"$U$", color="k", size=15)
    # x-axis
    pylab.arrow(0, 0, 1.05*max(r), 0, 
                fc='k', head_width=0.12, head_length=0.02)

    pylab.text(1.05*max(r), -0.025*max(U), r"$r$", color="k", size=15, verticalalignment="top")

    # draw a line representing a classical particle's energy
    E = 0.3*max(U)
    pylab.plot([0, max(r)], [E, E], 'b--')
    pylab.text(0.75*max(r), 1.05*E, "incoming proton KE", color='b',
               horizontalalignment="center")
    


    # find the x-coord of the classical turning point
    x = 1.0/E
    arrow_params = {'length_includes_head':True}

    pylab.arrow(x, 0, 0, E, fc='0.5', ec='0.5', head_width=0.015, head_length=0.25,
                **arrow_params)

    pylab.text(x, -0.025*max(U), r"$r_0$", color="0.5", size=15,horizontalalignment="center", verticalalignment="top")

    
    f = pylab.gcf()
    f.set_size_inches(6.0,6.0)

    pylab.savefig("coulomb_barrier.png")



if __name__== "__main__":
    coulomb_barrier()

