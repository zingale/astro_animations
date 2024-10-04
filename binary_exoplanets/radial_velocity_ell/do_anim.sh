#!/bin/bash

# make the movie
./radial_velocity_ell.py
ffmpeg -framerate 15 -f image2 -pattern_type glob -i "radial_velocity_????.png" -vcodec mpeg4 -c:v libx264 -crf 20 -pix_fmt yuv420p -movflags +faststart radial_velocity_ell.mp4

rm -f radial_velocity_????.png

