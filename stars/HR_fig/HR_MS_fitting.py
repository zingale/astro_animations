# plot a H-R diagram for a cluster for MS fitting.
#
# stellar properties approximate


import math
import numpy as np
import matplotlib.pyplot as plt

def HR_radius():

    # B-V plotting range
    BVmin = -0.5
    BVmax = 2.5

    # M plotting range
    Mmin = 20
    Mmax = -5


    #-------------------------------------------------------------------------
    # main sequence data 
    nstars = 13
    M  = np.zeros(nstars, np.float64)
    BV = np.zeros(nstars, np.float64)

    # B-V
    BV[:] = [-0.1, 0.0, 0.2, 0.3, 0.5, 0.8, 1.0, 1.2, 1.3, 1.4, 1.5, 1.6, 2.0]

    # absolute magnitude
    M[:] =  [0.5, 1,   2,   2.5, 4,   6,   6.75,7.5, 8,   8.75,11,  12,  16]

    #-------------------------------------------------------------------------
    ax = plt.subplot(111)
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    
    # draw and label the main sequence
    plt.plot(BV,M,'b-')     # true MS
    plt.plot(BV,M+5,'b:')   # cluster with unknown distance


    # reverse the x-axis
    plt.axis([BVmin, BVmax, Mmin, Mmax])
    ax = plt.gca()

    plt.xlabel("$B-V$")
    plt.ylabel("$M$")

    plt.subplots_adjust(left=0.125,right=0.95,bottom=0.1,top=0.95)

    f = plt.gcf()
    f.set_size_inches(6.0,7.0)

    plt.savefig("HR_MS_fitting.pdf")

if __name__== "__main__":
    HR_radius()

