#!/bin/bash

# make the movie
./ellipsedraw.py
ffmpeg -framerate 15 -f image2 -pattern_type glob -i "ellipsedraw_*.png" -vcodec mpeg4 -c:v libx264 -crf 20 -pix_fmt yuv420p -movflags +faststart ellipsedraw.mp4


# clean-up
rm -f ellipsedraw_???.png

