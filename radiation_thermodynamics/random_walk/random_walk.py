#!/bin/env python

import math
import random

import numpy as np
import pylab


class RandomWalk:

    def __init__(self, l, R, seed=1000):

        self.l = l
        self.R = R

        self.seed = seed

        self.max_steps = 5000

        self.x = None
        self.y = None


    def walk(self):

        # set the initial position to be the origin
        x_0 = 0.0
        y_0 = 0.0

        x = []
        y = []

        random.seed(self.seed)

        for n in range(self.max_steps):

            # compute a random angle
            angle = 2.0*math.pi*random.random()

            # compute the end coordinates of the segment
            x_1 = self.l*math.cos(angle) + x_0
            y_1 = self.l*math.sin(angle) + y_0

            x.append(x_1)
            y.append(y_1)
            
            # have we hit the edge of our domain?
            if math.sqrt(x_1**2 + y_1**2) >= self.R:
                break
            else:
                x_0 = x_1
                y_0 = y_1

        self.x = np.array(x)
        self.y = np.array(y)
            

def random_walk():

    # length of single step
    l = 1.0

    # maximum radius of domain
    R = 25.0

    # take steps, draw a segment in a random direction, and save the frame
    pylab.clf()

    # draw a circle to indicate the extent of the domain
    npts = 360
    theta = np.arange(npts)*2*math.pi/(npts-1)

    pylab.plot(R*np.cos(theta), R*np.sin(theta), color='k')

    pylab.subplots_adjust(left=0,right=1.0,bottom=0,top=1.0)

    # do a random walk
    r1 = RandomWalk(l, R, seed=1000)
    r1.walk()

    for n in range(len(r1.x)-1):
                
        pylab.plot([r1.x[n], r1.x[n+1]], 
                   [r1.y[n], r1.y[n+1]], color='r')

        pylab.axis([-1.1*R,1.1*R,-1.1*R,1.1*R])
        pylab.axis("off")

        f = pylab.gcf()
        f.set_size_inches(7.2, 7.2)

        pylab.savefig("random_walk_%04d.png" % n)


    # now do a second one
    r2 = RandomWalk(l, R, seed=1001)
    r2.walk()

    pylab.clf()

    pylab.plot(R*np.cos(theta), R*np.sin(theta), color='k')

    pylab.subplots_adjust(left=0,right=1.0,bottom=0,top=1.0)

    # draw the old one lighter
    for n in range(len(r1.x)-1):

        pylab.plot([r1.x[n], r1.x[n+1]], 
                   [r1.y[n], r1.y[n+1]], color='r', alpha=0.33)
        

    # now draw the new one
    for n in range(len(r2.x)-1):
                
        pylab.plot([r2.x[n], r2.x[n+1]], 
                   [r2.y[n], r2.y[n+1]], color='b')

        pylab.axis([-1.1*R,1.1*R,-1.1*R,1.1*R])
        pylab.axis("off")

        f = pylab.gcf()
        f.set_size_inches(7.2, 7.2)

        pylab.savefig("random_walk_%04d.png" % (n + len(r1.x)) )


if __name__== "__main__":
    random_walk()

