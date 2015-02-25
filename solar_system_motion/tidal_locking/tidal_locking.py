import matplotlib.pyplot as plt
import numpy as np

import anim_solvers.stick_figure as sf


class Circle:
    def __init__(self, R):
        theta = np.radians(np.arange(361))
        self.R = R
        self.x = R*np.cos(theta)
        self.y = R*np.sin(theta)
        


def doit():


    N = 360

    P_orbit = 1.0
    P_rotate = 0.1

    R_orbit = 25.0


    orbit = Circle(R_orbit)
    
    omega_orbit = 2.0*np.pi/P_orbit
    omega_rotate = 2.0*np.pi/P_rotate

    R_moon = 2.5
    R_planet = 10.0

    L = 1.5   # height of person
    
    moon = Circle(R_moon)
    planet = Circle(R_planet)    

    dt = P_orbit/float(N)
    
    for n in range(N):

        t = n*dt
        
        x_m = R_orbit*np.cos(omega_orbit*t)
        y_m = R_orbit*np.sin(omega_orbit*t)

        plt.clf()

        # draw the big planet
        plt.fill(planet.x, planet.y, color="b")

        # draw the orbit
        plt.plot(orbit.x, orbit.y, color="0.5", ls=":")

        # draw the moon
        plt.fill(x_m + moon.x, y_m + moon.y, color="0.5")

        # draw the person
        center = ( x_m + (R_moon + 0.5*L)*np.cos(omega_rotate*t),
                   y_m + (R_moon + 0.5*L)*np.sin(omega_rotate*t) )
        sf.draw_person(center, L, omega_rotate*t-np.pi/2, color="r")
        
        plt.axis("off")

        plt.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

        ax = plt.gca()
        ax.set_aspect("equal", "datalim")
        
        f = plt.gcf()
        f.set_size_inches(7.2,7.2)
                                                         
        plt.savefig("tidal_locking_{:04}.png".format(n))
        


if __name__ == "__main__":
    doit()
