#!/bin/bash

# make the movie
./planetary_orbits.py
ffmpeg -framerate 15 -f image2 -pattern_type glob -i "*.png" -vcodec mpeg4 -c:v libx264 -crf 20 -pix_fmt yuv420p -movflags +faststart planetary_orbits.mp4

rm -f planetary_orbit_????.png

