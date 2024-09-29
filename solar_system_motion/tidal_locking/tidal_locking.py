import matplotlib.pyplot as plt
import numpy as np

import anim_solvers.stick_figure as sf


class Circle:
    def __init__(self, R):
        theta = np.radians(np.arange(361))
        self.R = R
        self.x = R*np.cos(theta)
        self.y = R*np.sin(theta)
        

class Ellipse:
    def __init__(self, a, e, rot, theta_range=None):
        if theta_range == None:
            theta = np.radians(np.arange(361))
        else:
            # should probably use linspace...
            theta = np.radians(np.arange(theta_range[0],theta_range[1]+1))
            
        self.x = np.zeros_like(theta)
        self.y = np.zeros_like(theta)
        for n in range(len(theta)):
            self.x[n], self.y[n] = ellipse_point(a, e, theta[n], rot)
        

def ellipse_point(a, e, theta, rot):
    """ return the x, y coordinates of a point on a ellipse an angle
        theta from the x-axis, where theta is measured from the *center*
        of the ellipse.  Then rotate the point about the center by
        an angle rot"""

    # semi-minor axis
    b = a*np.sqrt(1 - e**2)

    # polar form relative to center
    r = a*b/np.sqrt( (b*np.cos(theta))**2 + (a*np.sin(theta))**2 )

    x = r*np.cos(theta)
    y = r*np.sin(theta)

    return sf._rotate([x, y], (0.,0.), rot)

            
def doit():


    N = 720

    P_orbit = 1.0
    P_rotate = 0.2

    R_orbit = 25.0


    orbit = Circle(R_orbit)
    
    omega_orbit = 2.0*np.pi/P_orbit

    R_moon = 2.5
    R_planet = 10.0

    L = 2.0   # height of person
    
    moon = Circle(R_moon)

    planet = Circle(R_planet)    

        
    dt = P_orbit/float(N)

    for d, scene in enumerate(["unlocked", "locked"]):

        if scene == "unlocked":
            omega_rotate = 2.0*np.pi/P_rotate
        else:
            omega_rotate = 2.0*np.pi/P_orbit
            
        
        for n in range(N):

            t = n*dt
            
            x_m = R_orbit*np.cos(omega_orbit*t)
            y_m = R_orbit*np.sin(omega_orbit*t)

            plt.clf()

            # draw the big planet
            plt.fill(planet.x, planet.y, color="#e6e6b8")

            # draw the orbit
            plt.plot(orbit.x, orbit.y, color="0.5", ls=":")

            # draw the moon
            #plt.fill(x_m + moon.x, y_m + moon.y, color="0.5")

            # draw the tides
            a = 1.4*R_moon
            e = 0.6

            # color the ellipse in segments to show the rotation more clearly
            n_segments = 12 
            for l, q in enumerate(zip(np.arange(n_segments)*360.0/n_segments, (np.arange(n_segments)+1)*360.0/n_segments)):
                theta_range = (q[0] + np.degrees((omega_rotate-omega_orbit)*t),
                               q[1] + np.degrees((omega_rotate-omega_orbit)*t))
                tides = Ellipse(a, e, omega_orbit*t, theta_range=theta_range)
                if l % 2 == 0:
                    c = "m"
                else:
                    c = "0.5"
                
                plt.fill([x_m] + list(x_m + tides.x) + [x_m],
                         [y_m] + list(y_m + tides.y) + [y_m],
                         color=c, zorder=-100, alpha=0.5)
        
        
            # draw the person

            # find the point on the elliptical moon where the person
            # lies.  The ellipse is already rotated by omega_orbit*t
            # and the person will be omega_rotate*t wrt the x-axis, so
            # they are (omega_rotate - omega_orbit)*t from the
            # semi-major axis of the ellipse.  So we use the polar
            # equation of an ellipse (wrt center) to find that point
            # and then rotate that point about the ellipse center to
            # correspond to the current orientation of the ellipse
            xp, yp = ellipse_point(a, e,
                                   (omega_rotate-omega_orbit)*t, omega_orbit*t)

            # now we extend by the half-height of the person along the
            # line connecting the person to the center of the ellipse,
            # since we draw a person by specifying the middle of the
            # person.
            xp1 = xp + L/2*np.cos(omega_rotate*t)
            yp1 = yp + L/2*np.sin(omega_rotate*t)
        
            center = (x_m + xp1, y_m + yp1)        
            sf.draw_person(center, L, omega_rotate*t-np.pi/2, color="r")

            # and a line to the Sun and Moon
            plt.plot([0, x_m], [0, y_m], "k", ls="--")
        
            plt.axis("off")

            plt.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95)

            max_dist = R_orbit + R_moon + 1.5*L
            plt.axis([-max_dist, max_dist, -max_dist, max_dist])

            f = plt.gcf()
            
            if scene == "unlocked":
                plt.text(0.5, 0.95,
                         "no tidal locking",
                         horizontalalignment="center",
                         transform=f.transFigure)
                plt.text(0.5, 0.05,
                         "the moon's rotation is independent of its orbital period around the planet",
                         horizontalalignment="center",
                         transform=f.transFigure)
                
            else:
                plt.text(0.5, 0.95,                
                         "tidally locked moon",
                         horizontalalignment="center",
                         transform=f.transFigure)
                plt.text(0.5, 0.05,
                         "the moon's rotation period is equal to its orbital period",
                         horizontalalignment="center",
                         transform=f.transFigure)
                
            ax = plt.gca()
            ax.set_aspect("equal", "datalim")
        

            f.set_size_inches(7.2,7.2)
    
            plt.savefig(f"tidal_locking_{n+d*N:04}.png")
        


if __name__ == "__main__":
    doit()
