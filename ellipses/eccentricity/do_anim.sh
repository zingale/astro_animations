#!/bin/bash

# make the movie
./eccentricity.py
mkmovie.py -o eccentricity ellipse_???.png

# store one of the frames as the preview frame for the webpage
cp ellipse_050.png eccentricity_preview.png

# clean-up
rm -f ellipse_???.png

