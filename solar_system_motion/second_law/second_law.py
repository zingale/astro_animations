#!/bin/env python

import math
import numpy
import pylab

import anim_solvers.solar_system_integrator as solar_system_integrator

# illustrate Kepler's second law by showing the areas corresponding to
# equal time intervals.

# M. Zingale 

# we work in units of AU, solar masses, and years

class interval:

    def __init__ (self, start, end):

        self.start = start
        self.end = end


def doit():

    s = solar_system_integrator.SolarSystem()

    # set the semi-major axis and eccentricity
    a = 1.5874
    e = 0.5

    # set the initial coordinates -- perihelion
    s.add_planet(a, e, loc="perihelion")

    # compute the period of the orbit from Kepler's law and make
    # the timestep by 1/720th of a period
    P = math.sqrt(a**3)

    # put 720 steps per period (the input is steps per year)
    sol = s.integrate(720/P, P)

    
    # set up some intervals (start time, end time) to shade.  
    # Note, the length of the intervals should be the same
    intervals = []
    intervals.append(interval(0, P/12))
    intervals.append(interval(P/2, 7.0*P/12.0))
    intervals.append(interval(9.5*P/12, 10.5*P/12.0))

    for i in range(len(intervals)):

        if intervals[i].start < 0 or intervals[i].end > P:
            print "ERROR: interval {} not contained in a single orbit".format(i)

        print "interval {}, dt = {}".format(i, intervals[i].end - intervals[i].start)


    # plot the orbit
    iframe = 0

    for n in range(len(sol[0].x)):

        pylab.clf()

        pylab.subplots_adjust(left=0.1,right=0.9,bottom=0.1,top=0.9)

        ax = pylab.gca()
        ax.set_aspect("equal", "datalim")
        pylab.axis("off")

        # plot the foci
        pylab.scatter([0], [0], s=350, marker=(5,1), color="k", zorder=10)
        pylab.scatter([0], [0], s=300, marker=(5,1), color="y", zorder=10)

        # plot planet 
        pylab.plot(sol[0].x, sol[0].y, color="r")
        pylab.scatter([sol[0].x[n]], [sol[0].y[n]], s=100, color="r", zorder=10)
        
        # shade the intervals
        for i in range(len(intervals)):

            # one vertex is the Sun
            x = [0]
            y = [0]

            for m in range(len(sol[0].x)):

                if sol[0].t[m] >= intervals[i].start and \
                   sol[0].t[m] <= min(intervals[i].end,sol[0].t[n]):
                    x.append(sol[0].x[m])
                    y.append(sol[0].y[m])

                elif (sol[0].t[m] > intervals[i].end):
                    break


            pylab.fill(x, y, alpha=0.20, facecolor="b")


        xmin = -a*(1.0+e)
        xmax = a*(1.0-e)
        dx = xmax-xmin
        pylab.axis([xmin-0.1*dx, xmax+0.1*dx, -1.5, 1.5])        

        f = pylab.gcf()
        f.set_size_inches(9.6,7.2)

        pylab.savefig("second_law_%04d.png" % n)

    
if __name__== "__main__":
    doit()


    
        
