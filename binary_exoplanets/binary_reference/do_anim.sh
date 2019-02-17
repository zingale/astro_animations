#!/bin/bash

# make the movie
./binary_reference.py
mkmovie.py -o binary_reference_mratio=4_e=0.4 binary_star_????.png

# store one of the frames as the preview frame for the webpage and clean-up
cp binary_star_0150.png binary_reference_mratio=4_e=0.4_preview.png
rm -f binary_star_????.png

