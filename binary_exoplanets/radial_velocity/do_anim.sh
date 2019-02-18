#!/bin/bash

# make the movie
./radial_velocity.py
mkmovie.py -o radial_velocity radial_velocity_????.png

# store one of the frames as the preview frame for the webpage and clean-up
cp radial_velocity_0450.png radial_velocity_preview.png
rm -f radial_velocity_????.png

