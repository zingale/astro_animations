import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy
import matplotlib as mpl


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
data[:, :] = 0.99
data[lonv > 180] = 0.1

# loop over the longitudes
i = 0
while i < len(lon_center):

    fig = plt.figure(figsize=(12.8, 7.2), facecolor="black")
    fig.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=0.95)


    # draw the sphere showing the phase
    ax = fig.add_axes([0.05, 0.333, 0.3, 0.333])  #], axisbg="k")

    bmap = Basemap(projection='ortho', lat_0=0, lon_0=lon_center[i] - 90,
                   resolution='l', area_thresh=1000., ax=ax)

    bmap.drawmapboundary()
    bmap.pcolor(lonv, latv, data, latlon=True, cmap=mpl.colormaps["bone"],
                vmin=0.0, vmax=1.0)

    ax.text(0.5, -0.2, "phase seen from Earth",
            horizontalalignment="center", color="w",
            transform=ax.transAxes)

    # draw the orbit seen top-down -- we adjust the phase angle here to sync
    # up with the phase view
    ax = fig.add_axes([0.4, 0.05, 0.55, 0.9])

    theta = numpy.linspace(0, 2.0*numpy.pi, 361, endpoint=True)
    theta_half = numpy.linspace(0, numpy.pi, 361, endpoint=True) - numpy.pi/2.0

    x = numpy.cos(theta)
    y = numpy.sin(theta)

    x_half = numpy.cos(theta_half)
    y_half = numpy.sin(theta_half)

    # draw the Earth
    # draw the Earth -- full circle (dark)
    ax.fill(0.125*x, 0.125*y, "0.25", zorder=100)

    # semi-circle -- illuminated
    ax.fill(0.125*x_half, 0.125*y_half, "C0", zorder=101)

    # draw the orbit
    ax.plot(x, y, "w--")

    # draw the Moon -- full circle (dark)
    ax.fill(numpy.cos(numpy.radians(lon_center[i])) + 0.075*x,
            numpy.sin(numpy.radians(lon_center[i])) + 0.075*y, "0.25",
            zorder=100)

    # semi-circle -- illuminated
    ax.fill(numpy.cos(numpy.radians(lon_center[i])) + 0.075*x_half,
            numpy.sin(numpy.radians(lon_center[i])) + 0.075*y_half, "1.0",
            zorder=101)

    # sunlight
    for ypos in [-0.6, -0.2, 0.2, 0.6]:
        ax.arrow(1.6, ypos, -0.4, 0.0, color="y",
                 length_includes_head=True,
                 head_width=0.1, width=0.05, overhang=-0.1)

    ax.text(1.7, 0.0, "sunlight", rotation=90, fontsize=16,
            horizontalalignment="center", verticalalignment="center",
            color="y")

    ax.axis([-1.1, 1.8, -1.1, 1.1])

    ax.set_aspect("equal", "datalim")

    ax.axis("off")

    fig.savefig(f"phase-{i:03d}.png", facecolor=fig.get_facecolor())
    plt.close(fig)

    i += 1
