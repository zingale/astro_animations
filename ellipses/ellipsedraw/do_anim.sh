#!/bin/bash

# make the movie
./ellipsedraw.py
mkmovie.py -o ellipsedraw ellipsedraw_???.png

# store one of the frames as the preview frame for the webpage
cp ellipsedraw_120.png ellipsedraw_preview.png

# clean-up
rm -f ellipsedraw_???.png

