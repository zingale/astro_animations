import pylab
from mpl_toolkits.basemap import Basemap
import numpy
import matplotlib.cm as cm


# to show the phase we will paint a hemisphere white and the other
# black.  We will then change the longitude that we are centered above
# through the full 360 degress, showing the range of phases.

# the range of longitudes to view
lon_center = numpy.linspace(0, 360, 361, endpoint=True)

# our grid of lon, lat that we will cover (half white, half dark gray)
# to show the illumination of the moon.
lon = numpy.linspace(0, 360, 360)
lat = numpy.linspace(-90, 90, 180)

lonv, latv = numpy.meshgrid(lon, lat)

data = numpy.zeros_like(lonv)
data[:,:] = 0.99
data[lonv > 180] = 0.1

# loop over the longitudes
i = 0
while i < len(lon_center):

    fig = pylab.figure(figsize=(12.8, 7.2), facecolor="black")


    #-------------------------------------------------------------------------
    # draw the sphere showing the phase
    pylab.axes([0.05, 0.333, 0.3, 0.333], axisbg="k")

    map = Basemap(projection='ortho', lat_0 = 0, lon_0 = lon_center[i] - 90,
                  resolution = 'l', area_thresh = 1000.)

    map.drawmapboundary()
    map.pcolor(lonv, latv, data, latlon=True, cmap=cm.bone, vmin = 0.0, vmax=1.0)

    ax = pylab.gca()
    pylab.text(0.5, -0.2, "phase seen from Earth", 
               horizontalalignment="center", color="w", 
               transform=ax.transAxes)

    #-------------------------------------------------------------------------
    # draw the orbit seen top-down -- we adjust the phase angle here to sync
    # up with the phase view 
    pylab.axes([0.4,0.1, 0.5, 0.8])

    theta = numpy.linspace(0, 2.0*numpy.pi, 361, endpoint=True)
    theta_half = numpy.linspace(0, numpy.pi, 361, endpoint=True) - numpy.pi/2.0

    x = numpy.cos(theta)
    y = numpy.sin(theta)

    x_half = numpy.cos(theta_half)
    y_half = numpy.sin(theta_half)


    # draw the Earth
    pylab.scatter([0],[0],s=1600,marker="o",color="k")
    pylab.scatter([0],[0],s=1500,marker="o",color="b")

    # draw the orbit
    pylab.plot(x, y, "w--")

    # draw the Moon -- full circle (dark)
    pylab.fill(numpy.cos(numpy.radians(lon_center[i])) + 0.075*x, 
               numpy.sin(numpy.radians(lon_center[i])) + 0.075*y, "0.25", zorder=100)

    # semi-circle -- illuminated
    pylab.fill(numpy.cos(numpy.radians(lon_center[i])) + 0.075*x_half, 
               numpy.sin(numpy.radians(lon_center[i])) + 0.075*y_half, "1.0", zorder=101)

    
    # sunlight
    pylab.arrow(1.6, -0.6, -0.4, 0.0, color="y",
                length_includes_head=True,
                head_width = 0.1, width=0.05, overhang=-0.1)

    pylab.arrow(1.6, -0.2, -0.4, 0.0, color="y",
                length_includes_head=True,
                head_width = 0.1, width=0.05, overhang=-0.1)

    pylab.arrow(1.6, 0.2, -0.4, 0.0, color="y",
                length_includes_head=True,
                head_width = 0.1, width=0.05, overhang=-0.1)

    pylab.arrow(1.6, 0.6, -0.4, 0.0, color="y",
                length_includes_head=True,
                head_width = 0.1, width=0.05, overhang=-0.1)

    pylab.text(1.7, 0.0, "sunlight", rotation=90, fontsize=16,
               horizontalalignment="center", verticalalignment="center", color="y")

    ax = pylab.gca()
    ax.set_aspect("equal", "datalim")

    pylab.axis("off")


    pylab.savefig("phase-{:03d}.png".format(i), facecolor=fig.get_facecolor())
    pylab.close()

    i += 1
