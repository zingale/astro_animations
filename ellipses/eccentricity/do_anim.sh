#!/bin/bash

# make the movie
./eccentricity.py
ffmpeg -framerate 15 -f image2 -pattern_type glob -i "ellipse_*.png" -vcodec mpeg4 -c:v libx264 -crf 20 eccentricity.mp4

# clean-up
rm -f ellipse_???.png

