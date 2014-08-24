import math
import numpy as np
import matplotlib.pyplot as plt


# Demonstrate the principle of Parallax.  We will take Earth's orbit
# to be circular.
#
# We work in units of AU, yr, and M_sun in these units, G = 4 pi^2

# M. Zingale 



class ParallaxScene:
    """
    We'll treat the entire collection of the Earth/Sun, foreground
    star, and background star as an object.  The only real thing that we
    need to change from frame to frame is the location of Earth
    """

    def __init__(self):

        # start Earth on the x-axis, on the opposite side of the field of
        # stars we will reference -- we accomplish this through a phase
        self.phi = math.pi
        
        # number of steps per year (make this a number divisible by 4)
        self.nsteps_year = 360   

        # angular velocity (radians per year)
        self.omega = 2.0*math.pi/1.0          

        # semi-major axis of planet Earth
        self.a_E = 1.0       

        # position of Earth over the year
        omega_t = np.arange(self.nsteps_year)*2.0*math.pi/(self.nsteps_year-1)
        self.x_orbit = self.a_E*np.cos(omega_t + self.phi)
        self.y_orbit = self.a_E*np.sin(omega_t + self.phi)

        # foreground star
        self.x_fg = 3.5
        self.y_fg = 0.0


    def draw_sun_and_orbit(self):

        # draw the Sun
        plt.scatter([0], [0], s=250, marker=(5,1), color="k")
        plt.scatter([0], [0], s=200, marker=(5,1), color="y")

        # plot the orbit
        plt.plot(self.x_orbit, self.y_orbit, "b--")


    def draw_earth(self, time, connect_to_fg=0):

        x_E = self.a_E*math.cos(self.omega*time + self.phi)
        y_E = self.a_E*math.sin(self.omega*time + self.phi)

        # plot Earth
        plt.scatter([x_E], [y_E], s=100, color="b")

        # draw the line connecting Earth and the foreground star
        if connect_to_fg == 1:
            slope = (y_E - self.y_fg)/(x_E - self.x_fg)
            xpt1 = 4.5
            ypt1 = y_E + slope*(xpt1 - x_E)
            x_E_old = x_E
            y_E_old = y_E
            plt.plot([x_E_old,xpt1], [y_E_old,ypt1], 'g--')


    def draw_foreground_star(self):
        # draw the foreground star
        plt.scatter([self.x_fg], [self.y_fg], s=200, marker=(5,1), color="r")        

        # draw the line connecting the Sun and the foreground star
        plt.plot([0,self.x_fg], [0,self.y_fg], 'k--')
        plt.text(1.5, -0.25, "d", color="k")


    def draw_background_stars(self):
        # draw some random background stars        

        pos = [(4.2, 2.1), (4.7, 1.0), (4.4, -0.4), (4.8, -0.9), 
               (4.1,-1.3), (4.3,-1.8), (4.5, 0.5)]

        for x, y in pos:
            plt.scatter( [x], [y], s=200, marker=(5,1), color="c")        



def parallax():


    t = 0.0

    nyears = 1.0

    # set the initial timestep
    nsteps_year = 360
    dt = 1.0/nsteps_year

    # compute the total number of steps needed 
    nsteps = int(nyears*nsteps_year)


    p = ParallaxScene()

    iout = 0

    # integrate until the Earth is at a right angle
    n = 0  
    while n < nsteps/4:

        plt.clf()
        
        p.draw_sun_and_orbit()
        p.draw_earth(t)
        p.draw_foreground_star()
        p.draw_background_stars()
        
        plt.axis([-1.5,5.0,-2.5,2.5])

        plt.axis('off')

        f = plt.gcf()
        f.set_size_inches(6.5,5.0)

        plt.xlabel("AU")
        plt.ylabel("AU")
        plt.text(-1.4,-2.0, "time = %6.3f yr" % t)
        plt.title("Parallax")
        
        plt.savefig("parallax_%04d.png" % iout)

        t += dt
        n += 1
        iout += 1




    # show the line connecting the current position and the foreground star
    # don't advance time
    print "connecting"

    nframes = 50

    plt.clf()

    p.draw_sun_and_orbit()
    p.draw_earth(t, connect_to_fg=1)
    p.draw_foreground_star()
    p.draw_background_stars()
    
    plt.axis([-1.5,5.0,-2.5,2.5])

    plt.axis('off')
        
    f = plt.gcf()
    f.set_size_inches(6.5,5.0)

    plt.xlabel("AU")
    plt.ylabel("AU")
    plt.text(-1.4,-2.0, "time = %6.3f yr" % t)
    plt.title("Parallax")

    plt.text(1.5,-0.8, "line of sight", color="g")
    plt.text(1.5,-1.0, "to foreground star", color="g")

    for n in range(nframes):
        plt.savefig("parallax_%04d.png" % iout)
        iout += 1


    # integrate for 6 months
    n = 0  
    while (n < nsteps/2):

        plt.clf()

        p.draw_sun_and_orbit()
        p.draw_earth(t)
        p.draw_foreground_star()
        p.draw_background_stars()

        plt.axis([-1.5,5.0,-2.5,2.5])

        plt.axis('off')

        f = plt.gcf()
        f.set_size_inches(6.5,5.0)

        plt.xlabel("AU")
        plt.ylabel("AU")
        plt.text(-1.4,-2.0, "time = %6.3f yr" % t)
        plt.title("Parallax")
        
        plt.savefig("parallax_%04d.png" % iout)

        t += dt
        n += 1
        iout += 1



    # show the new line connecting the current position and the foreground star
    # don't advance time
    print "connecting2"

    nframes = 50

    plt.clf()

    p.draw_sun_and_orbit()
    p.draw_earth(t, connect_to_fg=1)
    p.draw_foreground_star()
    p.draw_background_stars()

    plt.axis([-1.5,5.0,-2.5,2.5])

    plt.axis('off')

    f = plt.gcf()
    f.set_size_inches(6.5,5.0)

    plt.xlabel("AU")
    plt.ylabel("AU")
    plt.text(-1.4,-2.0, "time = %6.3f yr" % t)
    plt.title("Parallax")

    plt.text(1.5,1.0, "new line of sight", color="g")
    plt.text(1.5,0.8, "to foreground star", color="g")

    for n in range(nframes):
        plt.savefig("parallax_%04d.png" % iout)
        iout += 1


    # integrate for the final 1/4 year
    n = 0  
    while (n < nsteps/4):

        plt.clf()

        p.draw_sun_and_orbit()
        p.draw_earth(t)
        p.draw_foreground_star()
        p.draw_background_stars()

        plt.axis([-1.5,5.0,-2.5,2.5])

        plt.axis('off')

        f = plt.gcf()
        f.set_size_inches(6.5,5.0)

        plt.xlabel("AU")
        plt.ylabel("AU")
        plt.text(-1.4,-2.0, "time = %6.3f yr" % t)
        plt.title("Parallax")

        plt.savefig("parallax_%04d.png" % iout)

        t += dt
        n += 1
        iout += 1



    # summarize
    nframes = nsteps_year/2
    for n in range(nframes):

        plt.clf()
    
        p.draw_sun_and_orbit()
        p.draw_earth(t)
        p.draw_foreground_star()
        p.draw_background_stars()
        
        f = plt.gcf()
        f.set_size_inches(6.5,5.0)

        plt.xlabel("AU")
        plt.ylabel("AU")
        plt.text(-1.4,-2.0, "time = %6.3f yr" % t)
        plt.title("Parallax")
        
        plt.plot([0.0,0.0],[0.0,-1.0], "r")

        plt.text(2.5,-0.25,"p",color="g")
        plt.text(2.0,1.5,"tan p = 1 AU / d",color="g")
        plt.text(-0.5,-0.5, "1 AU", color="r")
        
        plt.axis([-1.5,5.0,-2.5,2.5])

        plt.axis('off')

        plt.savefig("parallax_%04d.png" % iout)

        t += dt
        iout += 1

    
if __name__== "__main__":
    parallax()


    
        
