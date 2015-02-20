#!/bin/env python

import math
import matplotlib.pyplot as plt
import random

import anim_solvers.solar_system_integrator as solar_system_integrator

# compute the orbit of Earth and Mars around the Sun.  Include a line
# connecting the two to show the line of sight, and therefore
# retrograde motion.

# M. Zingale


# extend our reference line (and star field) to both sides
both_sides = False


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


    # some background stars
    N = 10
    xpos = []
    ypos = []
    for s in range(N):
        xpos.append(random.uniform(3.0, 3.9))
        ypos.append(random.uniform(-1.4,1.9))

    if both_sides:
        xpos_left = []
        ypos_left = []
        for s in range(N):
            xpos_left.append(random.uniform(-2.4,-1.4))
            ypos_left.append(random.uniform(-1.4,1.9))


    # ================================================================
    # plotting
    # ================================================================
    for n in range(len(sol[0].x)):

        f = plt.figure()

        # plot the foci
        plt.scatter([0], [0], s=1600, marker=(20,1), color="k")
        plt.scatter([0], [0], s=1500, marker=(20,1), color="#FFFF00")

        # plot Earth
        plt.plot(sol[0].x, sol[0].y, color="b")
        plt.scatter([sol[0].x[n]], [sol[0].y[n]], s=100, color="b")

        # plot Mars
        plt.plot(sol[1].x, sol[1].y, color="r")
        plt.scatter([sol[1].x[n]],[sol[1].y[n]], s=100, color="r")

        # draw a line connecting Earth and Mars and extending a bit
        # further out
        slope = (sol[1].y[n] - sol[0].y[n])/(sol[1].x[n] - sol[0].x[n])
        xpt = 3.5
        ypt = sol[0].y[n] + slope*(xpt - sol[0].x[n])

        plt.plot([sol[0].x[n], xpt], [sol[0].y[n], ypt], "b--")

        # draw some random background stars
        for s in range(N):
            plt.scatter([xpos[s]], [ypos[s]], s=200, marker=(5,1), color="c")


        if both_sides:
            xpt = -2.5
            ypt = sol[0].y[n] + slope*(xpt - sol[0].x[n])
            plt.plot([sol[0].x[n], xpt], [sol[0].y[n], ypt], "b--")

            for s in range(N):
                plt.scatter([xpos_left[s]], [ypos_left[s]], s=200, marker=(5,1), color="c")


        plt.axis([-2.5, 4.0, -1.5, 2.0])
        plt.axis("off")

        ax = plt.gca()
        ax.set_aspect("equal", "datalim")

        f.set_size_inches(12.8,7.2)

        plt.xlabel("AU")
        plt.ylabel("AU")
        plt.text(-1.5,-1.5, "time = %6.3f yr" % sol[0].t[n])
        plt.title("Retrograde Mars")

        plt.savefig("retrograde_%04d.png" % n)

        plt.close()


if __name__== "__main__":
    doit()
