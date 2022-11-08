import math
import numpy as np
import matplotlib.pyplot as plt

import anim_solvers.solar_system_integrator as ssi

# compute the orbit of a 2 planets around the Sun using fourth-order accurate
# Runge-Kutta
#
# show two planets with the same perihelion, but different eccentricities
#
# M. Zingale (2010-09-16)

# we work in units of AU, yr, and M_sun
# in these units, G = 4 pi^2

# planet data
ecc_A = 0.0   # eccentricity of planet A
ecc_B = 0.4   # eccecntricity of planet B

a_A = 1.0     # semi-major axis of planet A

# we want the perihelion distance of planet B to be the same
# as planet a
r_p_A = a_A*(1.0 - ecc_A)
a_B = r_p_A/(1.0 - ecc_B)

ss = ssi.SolarSystem()
ss.add_planet(a_A, ecc_A, loc="perihelion")
ss.add_planet(a_B, ecc_B, loc="perihelion")

P_B = ss.period(1)



# integration data
nsteps_year = 365   # number of steps per year

sol = ss.integrate(nsteps_year, P_B)

sol_A = sol[0]
sol_B = sol[1]

nsteps = len(sol_A.x)

# plotting
plt.clf()

# plot the foci
plt.scatter([0],[0], s=350, marker=(20,1), color="k")
plt.scatter([0],[0], s=300, marker=(20,1), color="y")

# plot planet A orbit
plt.plot(sol_A.x, sol_A.y, color="r")

# plot planet B orbit
plt.plot(sol_B.x, sol_B.y, color="b", linestyle="--")

# plot planet B points
plt.scatter([sol_B.x[0]], [sol_B.y[0]], s=100, color="b", marker="h")
plt.text(sol_B.x[0]*1.15, sol_B.y[0], "B1", color="b")

plt.scatter([sol_B.x[nsteps//2]], [sol_B.y[nsteps//2]], s=100, color="b", marker="h")
plt.text(sol_B.x[nsteps//2]*1.15, sol_B.y[nsteps//2], "B2", color="b")

# the midpoint in time
x_mid = sol_B.x[3*nsteps//4]
y_mid = sol_B.y[3*nsteps//4]

plt.scatter([x_mid], [y_mid], s=100, color="b", marker="h")
plt.text(x_mid, y_mid*1.15, "B3", color="b")

# find the geometric 1/2 point between aphelion and perihelion -- this
# is just the semi-minor axis
b_B = a_B*np.sqrt(1.0 - ecc_B**2)

x_half = -a_B*ecc_B
y_half = -b_B

plt.scatter([x_half], [y_half], s=100, color="b", marker="h")
plt.text(x_half, y_half*1.15, "B4", color="b")

# last point equidistance spacing after the geometric 1/2 point
x_last = x_half + (x_half - x_mid)
y_last = y_mid

plt.scatter([x_last], [y_last], s=100, color="b", marker="h")
plt.text(x_last, y_last*1.15, "B5", color="b")    


# plot planet A points
plt.scatter([sol_A.x[0]],[sol_A.y[0]], s=100, color="r", marker="^")
plt.text(sol_A.x[0]*0.7, sol_A.y[0], "A1", color="r")



# draw the center of the ellipse
plt.scatter([-a_B*ecc_B], [0], s=100, color='b', marker='x')

# draw the semi-major and semi-minor axes
plt.plot([-a_B*ecc_B, -a_B*ecc_B], 
         [-a_B*math.sqrt(1.0 - ecc_B**2), a_B*math.sqrt(1.0 - ecc_B**2)], 
         color='b', linestyle=":")
plt.plot([-a_B*(1.0 + ecc_B), a_B*(1.0 - ecc_B)], [0,0], color='b', linestyle=":")

ax = plt.gca()
ax.set_aspect("equal", "datalim")
plt.axis("off")

#plt.axis([1.1*xmin, 1.5*xmax, -0.5*(xmax - xmin), 0.5*(xmax - xmin)])
plt.subplots_adjust(left=0.05, right=0.96, top=0.98, bottom=0.02)

f = plt.gcf()
f.set_size_inches(6.0,6.0)

plt.xlabel("AU")
plt.ylabel("AU")
    
outfile = "planets.pdf"
plt.savefig(outfile)



    
        
