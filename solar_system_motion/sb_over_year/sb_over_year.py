#!/bin/env python

# show the day night line over the course of a year (each frame is 2
# hours apart)
# centered on NY

from datetime import datetime, timedelta
from pytz import timezone
import numpy as np

from mpl_toolkits.basemap import Basemap

import matplotlib.pyplot as plt

#proj = "ortho"
proj = "moll"


# Stony Brook data 
lat_sb = 40.906
lon_sb = -73.128

if proj == "ortho":
    lat = lat_sb
    lon = lon_sb
else:
    lat = 0
    lon = 0


# work in UTC (+5)
date = datetime(2014, 01, 01, 12+5, 0, 0)
delta = timedelta(1.0)  # first argument is days

utc = timezone("UTC")
eastern = timezone('US/Eastern')


for n in range(365):

    plt.clf()

    map = Basemap(projection=proj, lat_0 = lat, lon_0 = lon,
                  resolution = 'l', area_thresh = 1000.)

    map.drawmapboundary()

    map.drawmeridians(np.arange(0, 360, 15), color="0.5", latmax=90)
    map.drawparallels(np.arange(-90, 90, 15), color="0.5", latmax=90)

    map.drawcoastlines()
    map.drawmapboundary(fill_color='aqua')
    map.fillcontinents(color='coral',lake_color='aqua')

    CS=map.nightshade(date)

    f = plt.gcf()
    f.set_size_inches(7.2, 7.2)

    plt.text(0.5, 0.95, "Noon at Stony Brook",
             transform=f.transFigure, horizontalalignment="center")

    plt.text(0.5, 0.05, "{}".format(utc.localize(date).astimezone(eastern)),
             transform=f.transFigure, horizontalalignment="center")

    plt.savefig("stonybrook_noon_{}_{:04}.png".format(proj, n))
    
    date += delta

