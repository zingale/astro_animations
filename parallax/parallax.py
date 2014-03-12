import math
import numpy
import pylab

# demonstrate the principle of Parallax.  We will take Earth's orbit to
# be circular.


# note: this routine is very klunky -- there is a lot of repitition that
# should be eliminated.

# M. Zingale (2008-09-02)

# we work in units of AU, yr, and M_sun
# in these units, G = 4 pi^2


# planet data
a_E = 1.0       # semi-major axis of planet Earth


# integration data
nsteps_year = 360   # number of steps per year 
                    # (make this a number divisible by 4)

nyears = 1.0        # number years to show

omega = 2.0*math.pi/1.0          # angular velocity (radians per year)

omega_orbit = numpy.arange(nsteps_year)*2.0*math.pi/(nsteps_year-1)
x_orbit = a_E*numpy.cos(omega_orbit)
y_orbit = a_E*numpy.sin(omega_orbit)

x_fg = 3.5
y_fg = 0.0


def parallax():


    t = 0.0
    dt = 0.0


    # ================================================================
    # Earth initialization
    # ================================================================

    # start Earth on the x-axis, on the opposite side of the field of
    # stars we will reference -- we accomplish this through a phase
    phi = math.pi


    # set the initial timestep
    dt = 1.0/nsteps_year

    # compute the total number of steps needed 
    nsteps = int(nyears*nsteps_year)

    print "nstep = ", nsteps

    iout = 0

    # integrate until the Earth is at a right angle
    n = 0  
    while (n < nsteps/4):

        x_E = a_E*math.cos(omega*t + phi)
        y_E = a_E*math.sin(omega*t + phi)

    
        # clear
        pylab.clf()


        # draw the line connecting the Sun and the foreground star
        pylab.plot([0,x_fg],[0,y_fg],'k--')
        pylab.text(1.5,-0.25,"d",color="k")

 
        # draw the Sun
        pylab.scatter([0],[0],s=250,marker=(5,1),color="k")
        pylab.scatter([0],[0],s=200,marker=(5,1),color="y")

        # draw the foreground star
        pylab.scatter([x_fg],[y_fg],s=200,marker=(5,1),color="r")        

        # plot the orbit
        pylab.plot(x_orbit, y_orbit, "b--")

        # plot Earth
        pylab.scatter([x_E],[y_E],s=100,color="b")

        # draw some random background stars
        pylab.scatter([4.2],[ 2.1],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.7],[ 1.0],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.4],[-0.4],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.8],[-0.9],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.1],[-1.3],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.3],[-1.8],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.5],[ 0.5],s=200,marker=(5,1),color="c")        

        

        pylab.axis([-1.5,5.0,-2.5,2.5])

        pylab.axis('off')

        f = pylab.gcf()
        f.set_size_inches(6.5,5.0)

        pylab.xlabel("AU")
        pylab.ylabel("AU")
        pylab.text(-1.4,-2.0, "time = %6.3f yr" % t)
        pylab.title("Parallax")


        outfile = "parallax_%04d.png" % iout
        pylab.savefig(outfile)

        t += dt
        n += 1
        iout += 1




    # show the line connecting the current position and the foreground star
    # don't advance time
    print "connecting"

    nframes = 50
    n = 0
    while (n < nframes):

        x_E = a_E*math.cos(omega*t + phi)
        y_E = a_E*math.sin(omega*t + phi)


        # clear
        pylab.clf()


        # draw the line connecting the Sun and the foreground star
        pylab.plot([0,x_fg],[0,y_fg],'k--')
        pylab.text(1.5,-0.25,"d",color="k")

 
        # draw the Sun
        pylab.scatter([0],[0],s=250,marker=(5,1),color="k")
        pylab.scatter([0],[0],s=200,marker=(5,1),color="y")

        # draw the line connecting Earth and the foreground star
        slope = (y_E - y_fg)/(x_E - x_fg)
        xpt1 = 4.5
        ypt1 = y_E + slope*(xpt1 - x_E)
        x_E_old = x_E
        y_E_old = y_E
        pylab.plot([x_E_old,xpt1],[y_E_old,ypt1], 'g--')

        # draw the foreground star
        pylab.scatter([x_fg],[y_fg],s=200,marker=(5,1),color="r")        

        # plot the orbit
        pylab.plot(x_orbit, y_orbit, "b--")

        # plot Earth        
        pylab.scatter([x_E],[y_E],s=100,color="b")

        # draw some random background stars
        pylab.scatter([4.2],[ 2.1],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.7],[ 1.0],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.4],[-0.4],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.8],[-0.9],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.1],[-1.3],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.3],[-1.8],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.5],[ 0.5],s=200,marker=(5,1),color="c")        

        

        pylab.axis([-1.5,5.0,-2.5,2.5])

        pylab.axis('off')

        f = pylab.gcf()
        f.set_size_inches(6.5,5.0)

        pylab.xlabel("AU")
        pylab.ylabel("AU")
        pylab.text(-1.4,-2.0, "time = %6.3f yr" % t)
        pylab.title("Parallax")

        pylab.text(1.5,-0.8, "line of sight", color="g")
        pylab.text(1.5,-1.0, "to foreground star", color="g")

        outfile = "parallax_%04d.png" % iout
        pylab.savefig(outfile)

        
        n += 1
        iout += 1




    # integrate for 6 months
    n = 0  
    while (n < nsteps/2):

        x_E = a_E*math.cos(omega*t + phi)
        y_E = a_E*math.sin(omega*t + phi)


        # clear
        pylab.clf()


        # draw the line connecting the Sun and the foreground star
        pylab.plot([0,x_fg],[0,y_fg],'k--')
        pylab.text(1.5,-0.25,"d",color="k")

 
        # draw the Sun
        pylab.scatter([0],[0],s=250,marker=(5,1),color="k")
        pylab.scatter([0],[0],s=200,marker=(5,1),color="y")

        # draw the line connecting old Earth and the foreground star
        pylab.plot([x_E_old,xpt1],[y_E_old,ypt1], 'g--')

        # draw the foreground star
        pylab.scatter([x_fg],[y_fg],s=200,marker=(5,1),color="r")        

        # plot the orbit
        pylab.plot(x_orbit, y_orbit, "b--")

        # plot Earth
        pylab.scatter([x_E],[y_E],s=100,color="b")

        # draw some random background stars
        pylab.scatter([4.2],[ 2.1],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.7],[ 1.0],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.4],[-0.4],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.8],[-0.9],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.1],[-1.3],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.3],[-1.8],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.5],[ 0.5],s=200,marker=(5,1),color="c")        

        

        pylab.axis([-1.5,5.0,-2.5,2.5])

        pylab.axis('off')

        f = pylab.gcf()
        f.set_size_inches(6.5,5.0)

        pylab.xlabel("AU")
        pylab.ylabel("AU")
        pylab.text(-1.4,-2.0, "time = %6.3f yr" % t)
        pylab.title("Parallax")


        outfile = "parallax_%04d.png" % iout
        pylab.savefig(outfile)

        t += dt
        n += 1
        iout += 1



    # show the new line connecting the current position and the foreground star
    # don't advance time
    print "connecting2"

    nframes = 50
    n = 0
    while (n < nframes):

        x_E = a_E*math.cos(omega*t + phi)
        y_E = a_E*math.sin(omega*t + phi)


        # clear
        pylab.clf()


        # draw the line connecting the Sun and the foreground star
        pylab.plot([0,x_fg],[0,y_fg],'k--')
        pylab.text(1.5,-0.25,"d",color="k")

 
        # draw the Sun
        pylab.scatter([0],[0],s=250,marker=(5,1),color="k")
        pylab.scatter([0],[0],s=200,marker=(5,1),color="y")

        # draw the line connecting Earth and the foreground star
        slope = (y_E - y_fg)/(x_E - x_fg)
        xpt2 = 4.5
        ypt2 = y_E + slope*(xpt2 - x_E)
        x_E_new = x_E
        y_E_new = y_E
        pylab.plot([x_E_new,xpt2],[y_E_new,ypt2], 'g--')
        pylab.plot([x_E_old,xpt1],[y_E_old,ypt1], 'g--')


        # draw the foreground star
        pylab.scatter([x_fg],[y_fg],s=200,marker=(5,1),color="r")        

        # plot the orbit
        pylab.plot(x_orbit, y_orbit, "b--")

        # plot Earth        
        pylab.scatter([x_E],[y_E],s=100,color="b")

        # draw some random background stars
        pylab.scatter([4.2],[ 2.1],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.7],[ 1.0],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.4],[-0.4],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.8],[-0.9],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.1],[-1.3],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.3],[-1.8],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.5],[ 0.5],s=200,marker=(5,1),color="c")        

        

        pylab.axis([-1.5,5.0,-2.5,2.5])

        pylab.axis('off')

        f = pylab.gcf()
        f.set_size_inches(6.5,5.0)

        pylab.xlabel("AU")
        pylab.ylabel("AU")
        pylab.text(-1.4,-2.0, "time = %6.3f yr" % t)
        pylab.title("Parallax")

        pylab.text(1.5,1.0, "new line of sight", color="g")
        pylab.text(1.5,0.8, "to foreground star", color="g")

        outfile = "parallax_%04d.png" % iout
        pylab.savefig(outfile)

        n += 1
        iout += 1



    # integrate for the final 1/4 year
    n = 0  
    while (n < nsteps/4):

        # clear
        pylab.clf()

        # draw the line connecting the Sun and the foreground star
        pylab.plot([0,x_fg],[0,y_fg],'k--')
        pylab.text(1.5,-0.25,"d",color="k")

 
        # draw the Sun
        pylab.scatter([0],[0],s=250,marker=(5,1),color="k")
        pylab.scatter([0],[0],s=200,marker=(5,1),color="y")

        # draw the foreground star
        pylab.scatter([x_fg],[y_fg],s=200,marker=(5,1),color="r")        

        # draw the connecting lines
        pylab.plot([x_E_new,xpt2],[y_E_new,ypt2], 'g--')
        pylab.plot([x_E_old,xpt1],[y_E_old,ypt1], 'g--')

        # plot the orbit
        pylab.plot(x_orbit, y_orbit, "b--")

        # plot Earth
        x_E = a_E*math.cos(omega*t + phi)
        y_E = a_E*math.sin(omega*t + phi)

        
        pylab.scatter([x_E],[y_E],s=100,color="b")

        # draw some random background stars
        pylab.scatter([4.2],[ 2.1],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.7],[ 1.0],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.4],[-0.4],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.8],[-0.9],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.1],[-1.3],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.3],[-1.8],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.5],[ 0.5],s=200,marker=(5,1),color="c")        

        

        pylab.axis([-1.5,5.0,-2.5,2.5])

        pylab.axis('off')

        f = pylab.gcf()
        f.set_size_inches(6.5,5.0)

        pylab.xlabel("AU")
        pylab.ylabel("AU")
        pylab.text(-1.4,-2.0, "time = %6.3f yr" % t)
        pylab.title("Parallax")


        outfile = "parallax_%04d.png" % iout
        pylab.savefig(outfile)

        t += dt
        n += 1
        iout += 1



    # summarize
    nframes = nsteps_year/2
    n = 0
    while (n < nframes):

        x_E = a_E*math.cos(omega*t + phi)
        y_E = a_E*math.sin(omega*t + phi)

        # clear
        pylab.clf()

        # draw the line connecting the Sun and the foreground star
        pylab.plot([0,x_fg],[0,y_fg],'k--')
        pylab.text(1.5,-0.25,"d",color="k")

 
        # draw the Sun
        pylab.scatter([0],[0],s=250,marker=(5,1),color="k")
        pylab.scatter([0],[0],s=200,marker=(5,1),color="y")

        # draw the foreground star
        pylab.scatter([x_fg],[y_fg],s=200,marker=(5,1),color="r")        

        # draw the connecting lines
        pylab.plot([x_E_new,xpt2],[y_E_new,ypt2], 'g--')
        pylab.plot([x_E_old,xpt1],[y_E_old,ypt1], 'g--')

        # plot the orbit
        pylab.plot(x_orbit, y_orbit, "b--")

        
        # plot Earth
        pylab.scatter([x_E],[y_E],s=100,color="b")

        # draw some random background stars
        pylab.scatter([4.2],[ 2.1],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.7],[ 1.0],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.4],[-0.4],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.8],[-0.9],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.1],[-1.3],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.3],[-1.8],s=200,marker=(5,1),color="c")        
        pylab.scatter([4.5],[ 0.5],s=200,marker=(5,1),color="c")        

        
        f = pylab.gcf()
        f.set_size_inches(6.5,5.0)

        pylab.xlabel("AU")
        pylab.ylabel("AU")
        pylab.text(-1.4,-2.0, "time = %6.3f yr" % t)
        pylab.title("Parallax")

        pylab.plot([0.0,0.0],[0.0,-1.0], "r")

        pylab.text(2.5,-0.25,"p",color="g")
        pylab.text(2.0,1.5,"tan p = 1 AU / d",color="g")
        pylab.text(-0.5,-0.5, "1 AU", color="r")

        pylab.axis([-1.5,5.0,-2.5,2.5])

        pylab.axis('off')

        outfile = "parallax_%04d.png" % iout
        pylab.savefig(outfile)

        t += dt
        n += 1
        iout += 1



    
if __name__== "__main__":
    parallax()


    
        
