import pylab
from mpl_toolkits.basemap import Basemap
import numpy
import math
import matplotlib.cm as cm


lon_center = numpy.linspace(0, 360, 25, endpoint=True)

print lon_center

lon = numpy.linspace(0, 360, 360)
lat = numpy.linspace(-90, 90, 180)

lonv, latv = numpy.meshgrid(lon, lat)

data = numpy.zeros_like(lonv)
data[lonv > 180] = 1.0

i = 0
while i < len(lon_center):

    map = Basemap(projection='ortho', lat_0 = 0, lon_0 = lon_center[i],
                  resolution = 'l', area_thresh = 1000.)

    map.drawmapboundary()

    #map.drawmeridians(numpy.arange(0, 360, 15), color="0.5", latmax=90)
    #map.drawparallels(numpy.arange(-90, 90, 15), color="0.5", latmax=90) 

    map.pcolor(lonv, latv, data, latlon=True, cmap=cm.binary)

    pylab.savefig("phase-{:02d}.png".format(i))

    i += 1
