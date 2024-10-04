#!/bin/env python

# show the day night line over the course of a year (each frame is 2
# hours apart)
# centered on NY

from datetime import datetime, timedelta
from pytz import timezone
import numpy as np

from mpl_toolkits.basemap import Basemap

import matplotlib.pyplot as plt

def stonybrook(proj="ortho"):

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
    date = datetime(2024, 1, 1, 12+5, 0, 0)
    delta = timedelta(1.0)  # first argument is days

    utc = timezone("UTC")
    eastern = timezone('US/Eastern')

    for n in range(365):

        fig, ax = plt.subplots()

        map = Basemap(projection=proj, lat_0=lat, lon_0=lon,
                      resolution='l', area_thresh= 1000., ax=ax)

        map.drawmapboundary()

        map.drawmeridians(np.arange(0, 360, 15), color="0.5", latmax=90)
        map.drawparallels(np.arange(-90, 90, 15), color="0.5", latmax=90)

        map.drawcoastlines()
        map.drawmapboundary(fill_color='aqua')
        map.fillcontinents(color='coral',lake_color='aqua')

        CS=map.nightshade(date)

        fig.set_size_inches(7.2, 7.2)

        ax.text(0.5, 0.95, "Noon at Stony Brook",
                transform=fig.transFigure, horizontalalignment="center", fontsize=20)

        ax.text(0.5, 0.05, "{}".format(utc.localize(date).astimezone(eastern)),
                transform=fig.transFigure, horizontalalignment="center", fontsize=16)

        fig.savefig(f"stonybrook_noon_{proj}_{n:04}.png")
        plt.close(fig)

        date += delta


if __name__ == "__main__":
    stonybrook(proj="ortho")
    stonybrook(proj="moll")
