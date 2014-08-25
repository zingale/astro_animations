#!/bin/env python

import math
import numpy
import pylab

import anim_solvers.solar_system_integrator as solar_system_integrator

# compute the orbit of Earth and Mars around the Sun.  Include a line
# connecting the two to show the line of sight, and therefore
# retrograde motion.

# M. Zingale

def doit():

    # planet data
    ecc_E = 0.016710  # eccentricity of planet Earth
    ecc_M = 0.093315  # eccecntricity of planet Mars

    a_E = 1.0       # semi-major axis of planet Earth
    a_M = 1.523679  # semi-major axis of planet Mars

    # integration data
    nsteps_year = 365*2   # number of steps per year
    nyears = 0.6          # total integration time (years)

    s = solar_system_integrator.SolarSystem()


    # Earth initialization

    # set the initial conditions.  The initial position is perihelion
    #y[0] = a_E*(1.0 - ecc_E)   # x 
    #y[1] = 0.0                 # y

    # at perihelion, all the veloicity is in the y-direction.
    #y[2] = 0.0      # v_x

    # v_y^2 = (GM/a) (1+e)/(1-e)  (see C&O Eq. 2.33 for example)
    # This is the perihelion velocity.  For a = 1, e = 0, v = 2 pi
    #y[3] = -math.sqrt( (G*M_sun/a_E) * (1.0 + ecc_E) / (1.0 - ecc_E))


    # Mars initialization

    # set the initial conditions.  The initial position is perihelion
    #y[4] = a_M*(1.0 - ecc_M)  # x 
    #y[5] = 0.0                 # y
    
    # at perihelion, all the veloicity is in the y-direction.
    #y[6] = 0.0      # v_x
    
    # v_y^2 = (GM/a) (1+e)/(1-e)  (see C&O Eq. 2.33 for example)
    # This is the perihelion velocity.  For a = 1, e = 0, v = 2 pi
    #y[7] = -math.sqrt( (G*M_sun/a_M) * (1.0 + ecc_M) / (1.0 - ecc_M))
    

    # These initial conditions were found by putting Mars and Earth
    # both at their perihelion, but with their velocities in the
    # opposite direction (i.e. we want them to go backwards).  This
    # configuration is opposition, and is setup by using the commented
    # out initial conditions above.  We then integrated backwards for 1/4
    # year (91 steps) to get these starting coordinates (note, we reverse
    # the direction of the velocity to get it orbitting in the forward
    # direction.)
    

    # Earth
    x  = -0.04631900088483218
    y  = -0.9994219951994862
    vx = 6.277324691390798
    vy = -0.185920887199495

    vp0 = solar_system_integrator.PlanetPosVel(x, y, vx, vy)
    s.add_planet(a_E, ecc_E, loc="specify", pos_vel=vp0)


    # Mars
    x  = 0.7856599524256417
    y  = -1.203323492875661
    vx = 4.280834571016523
    vy = 3.272064392180777

    vp0 = solar_system_integrator.PlanetPosVel(x, y, vx, vy)
    s.add_planet(a_M, ecc_M, loc="specify", pos_vel=vp0)
    
    
    # integrate
    sol = s.integrate(nsteps_year, nyears)


    # ================================================================
    # plotting
    # ================================================================
    for n in range(len(sol[0].x)):
        
        f = pylab.figure()

        # plot the foci
        pylab.scatter([0], [0], s=250,marker=(5,1), color="k")
        pylab.scatter([0], [0], s=200,marker=(5,1), color="y")

        # plot Earth
        pylab.plot(sol[0].x, sol[0].y, color="b")
        pylab.scatter([sol[0].x[n]], [sol[0].y[n]], s=100, color="b")

        # plot Mars
        pylab.plot(sol[1].x, sol[1].y, color="r")
        pylab.scatter([sol[1].x[n]],[sol[1].y[n]], s=100, color="r")

        # draw a line connecting Earth and Mars and extending a bit
        # further out
        slope = (sol[1].y[n] - sol[0].y[n])/(sol[1].x[n] - sol[0].x[n])
        xpt = 3.5
        ypt = sol[0].y[n] + slope*(xpt - sol[0].x[n])

        pylab.plot([sol[0].x[n], xpt], [sol[0].y[n], ypt], "b--")

        # draw some random background stars
        pylab.scatter([3.2],[ 1.9],s=200,marker=(5,1),color="c")        
        pylab.scatter([3.6],[ 1.0],s=200,marker=(5,1),color="c")        
        pylab.scatter([3.4],[-0.4],s=200,marker=(5,1),color="c")        
        pylab.scatter([3.8],[-0.9],s=200,marker=(5,1),color="c")        
        pylab.scatter([3.1],[-1.3],s=200,marker=(5,1),color="c")        

        pylab.axis([-2.,4.,-1.5,2.0])
        pylab.axis("off")

        ax = pylab.gca()
        ax.set_aspect("equal", "datalim")

        f.set_size_inches(12.8,7.2)

        pylab.xlabel("AU")
        pylab.ylabel("AU")
        pylab.text(-1.5,-1.5, "time = %6.3f yr" % sol[0].t[n])
        pylab.title("Retrograde Mars")

        pylab.savefig("retrograde_%04d.png" % n)

        pylab.close()

        n += 1

    
if __name__== "__main__":
    doit()


    
        
