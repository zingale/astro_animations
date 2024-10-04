#!/bin/bash

# make the movie
./binary_reference.py
ffmpeg -framerate 15 -f image2 -pattern_type glob -i "*.png" -vcodec mpeg4 -c:v libx264 -crf 20 binary_reference_mratio=4_e=0.4.mp4

rm -f binary_star_????.png

