#!/bin/bash

# make the movie
./radial_velocity.py
ffmpeg -framerate 15 -f image2 -pattern_type glob -i "radial_velocity_????.png" -vcodec mpeg4 -c:v libx264 -crf 20 radial_velocity.mp4

# store one of the frames as the preview frame for the webpage and clean-up
cp radial_velocity_0450.png radial_velocity_preview.png
rm -f radial_velocity_????.png

