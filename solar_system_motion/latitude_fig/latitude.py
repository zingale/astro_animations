#!/usr/bin/env python

import anim_solvers.stick_figure as sf
import numpy as np
import matplotlib.pylab as plt

def doit():

    # latitude
    lat = 40

    transparent = True
    
    # scalings for the person and earth
    L = 1.5
    R = 10

    theta = np.radians(np.arange(0,361))

    # draw a circle
    plt.plot(R*np.cos(theta), R*np.sin(theta), c="b")

    # draw equator and rotation axis
    plt.plot([-R, R], [0,0], c="k")
    plt.plot([0,0], [-R, R], c="k", ls="--")

    plt.text(-0.5*R, 0.25, "equator", horizontalalignment="center")
    
    # draw a person
    center = ( (R + 0.5*L)*np.cos(np.radians(lat)),
               (R + 0.5*L)*np.sin(np.radians(lat)) )
    sf.draw_person(center, L, np.radians(lat - 90), color="r")


    # draw the latitude angle line to zenith
    npts = 10
    x = np.arange(npts)*2*R/npts
    y = center[1]/center[0] * x

    plt.plot(x, y, "k", ls=":")

    plt.text(0.2*R, 0.5*0.2*R*np.sin(lat), r"$l$",
             fontsize=16, verticalalignment="center")
    
    plt.axis("off")
    
    ax = plt.gca()
    ax.set_aspect("equal", "datalim")
        
    plt.subplots_adjust(left=0.05, right=0.98, bottom=0.05, top=0.98)
    plt.axis([-1.2*R, 1.8*R, -1.2*R, 1.8*R])

    f = plt.gcf()
    f.set_size_inches(7.2, 7.2)
                            
    plt.savefig("latitude1.png", transparent=transparent)



    # add the zenith label
    plt.text(x[0.7*npts], y[0.7*npts], "points to zenith",
             horizontalalignment="center", verticalalignment="center",
             rotation = lat)

    plt.savefig("latitude2.png", transparent=transparent)


    # draw the local horizon -- it passed through center, perpendicular
    # to lat
    scenter = ( R*np.cos(np.radians(lat)),
                R*np.sin(np.radians(lat)) )

    p = (scenter[0] + R*np.cos(np.radians(90+lat)),
         scenter[1] + R*np.sin(np.radians(90+lat)) )

    yh = (p[1] - scenter[1])/(p[0] - scenter[0])*(x - scenter[0]) + scenter[1]

    plt.plot(x, yh, color="c")
    plt.text(x[0.7*npts], yh[0.7*npts], "local horizon",
             rotation=270+lat, color="c")             

    
    plt.savefig("latitude3.png", transparent=transparent)
    
    
    
if __name__ == "__main__":
    doit()

    
