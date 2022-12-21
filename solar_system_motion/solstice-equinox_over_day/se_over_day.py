#!/bin/env python

# show the day night line over the course of a year (each frame is 2
# hours apart)
# centered on NY

from datetime import datetime, timedelta
from pytz import timezone
import numpy as np

from mpl_toolkits.basemap import Basemap

import matplotlib.pyplot as plt

proj = "ortho"
#proj = "moll"


# Stony Brook data
lat_sb = 40.906
lon_sb = -73.128

if proj == "ortho":
    lat = lat_sb
    lon = lon_sb
else:
    lat = 0
    lon = 0


dmin = 2.5
N = int(24*60.0/dmin)

# work in UTC (+4 in summer for EDT)
date = datetime(2022, 12, 21, 0+4, 0, 0)
delta = timedelta(minutes=dmin)  # first argument is days

utc = timezone("UTC")
eastern = timezone('US/Eastern')


for n in range(N):

    plt.clf()

    bmap = Basemap(projection=proj, lat_0 = lat, lon_0 = lon,
                  resolution = 'l', area_thresh = 100.)

    bmap.drawmapboundary()

    bmap.drawmeridians(np.arange(0, 360, 15), color="0.5", latmax=90)
    bmap.drawparallels(np.arange(-90, 90, 15), color="0.5", latmax=90)

    bmap.drawparallels(np.array([-66.5, -23.5, 0, 23.5, 66.5]),
                      color="b", linewidth=2, latmax=90,
                      dashes=[1000,0.001])

    bmap.drawcoastlines()
    bmap.drawmapboundary(fill_color='aqua')
    bmap.fillcontinents(color='coral',lake_color='aqua')

    CS = bmap.nightshade(date)

    f = plt.gcf()
    f.set_size_inches(10.8, 10.8)

    plt.text(0.5, 0.95, "Stony Brook on the Winter Solstice",
             transform=f.transFigure, horizontalalignment="center")

    plt.text(0.5, 0.05, f"{utc.localize(date).astimezone(eastern)}",
             transform=f.transFigure, horizontalalignment="center")

    plt.savefig(f"sb_winter_solstice_{proj}_{n:04}.png")

    date += delta
