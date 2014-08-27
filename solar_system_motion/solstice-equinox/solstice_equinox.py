#!/bin/env python

from datetime import datetime
import numpy as np

from mpl_toolkits.basemap import Basemap

import matplotlib.pyplot as plt


ss_date = datetime(2014, 06, 21, 12, 0, 0)
eq_date = datetime(2014, 9, 23, 12, 0, 0)
ws_date = datetime(2014, 12, 21, 12, 0, 0)    

dates = [ss_date, eq_date, ws_date]
events = ["summer_solstice", "equinox", "winter_solstice"]
pretty_title = ["Summer Solstice", "Equinox", "Winter Solstice"]

for e, p, d in zip(events, pretty_title, dates):

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

        CS=map.nightshade(d)

        f = plt.gcf()
        f.set_size_inches(7.2, 7.2)

        plt.text(0.5, 0.95, p, transform=f.transFigure, 
                 horizontalalignment="center", fontsize="large")

        plt.text(0.5, 0.05, "noon UTC", transform=f.transFigure, 
                 horizontalalignment="center")

        plt.savefig("{}_{:03}.png".format(e, lon))

