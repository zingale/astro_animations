#!/bin/env python

from datetime import datetime
import numpy as np

from mpl_toolkits.basemap import Basemap

import matplotlib.pyplot as plt


# summer solstice
for lon in range(0, 360, 1):

    plt.clf()

    map = Basemap(projection='ortho', lat_0 = 0, lon_0 = lon,
                  resolution = 'l', area_thresh = 1000.)

    map.drawmapboundary()

    map.drawmeridians(np.arange(0, 360, 15), color="0.5", latmax=90)
    map.drawparallels(np.arange(-90, 90, 15), color="0.5", latmax=90)

    map.drawcoastlines()
    map.drawmapboundary(fill_color='aqua')
    map.fillcontinents(color='coral',lake_color='aqua')

    date = datetime(2014, 06, 21, 12, 0, 0)
    CS=map.nightshade(date)

    f = plt.gcf()
    f.set_size_inches(7.2, 7.2)

    plt.savefig("summer_solstice_{:03}.png".format(lon))


# equinox
for lon in range(0, 360, 1):

    plt.clf()

    map = Basemap(projection='ortho', lat_0 = 0, lon_0 = lon,
                  resolution = 'l', area_thresh = 1000.)

    map.drawmapboundary()

    map.drawmeridians(np.arange(0, 360, 15), color="0.5", latmax=90)
    map.drawparallels(np.arange(-90, 90, 15), color="0.5", latmax=90)

    map.drawcoastlines()
    map.drawmapboundary(fill_color='aqua')
    map.fillcontinents(color='coral',lake_color='aqua')

    date = datetime(2014, 9, 23, 12, 0, 0)
    CS=map.nightshade(date)

    f = plt.gcf()
    f.set_size_inches(7.2, 7.2)

    plt.savefig("equinox_{:03}.png".format(lon))


# winter solstice
for lon in range(0, 360, 1):

    plt.clf()

    map = Basemap(projection='ortho', lat_0 = 0, lon_0 = lon,
                  resolution = 'l', area_thresh = 1000.)

    map.drawmapboundary()

    map.drawmeridians(np.arange(0, 360, 15), color="0.5", latmax=90)
    map.drawparallels(np.arange(-90, 90, 15), color="0.5", latmax=90)

    map.drawcoastlines()
    map.drawmapboundary(fill_color='aqua')
    map.fillcontinents(color='coral',lake_color='aqua')

    date = datetime(2014, 12, 21, 12, 0, 0)
    CS=map.nightshade(date)

    f = plt.gcf()
    f.set_size_inches(7.2, 7.2)

    plt.savefig("winter_solstice_{:03}.png".format(lon))


