#!/bin/bash

# make the movie
./planetary_orbits.py
mkmovie.py -o planetary_orbits planetary_orbit_????.png

# store one of the frames as the preview frame for the webpage and clean-up
cp planetary_orbit_0150.png planetary_orbit_preview.png
rm -f planetary_orbit_????.png

