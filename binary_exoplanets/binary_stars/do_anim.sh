#!/bin/bash

# ---- M1/M2 = 1.0; e = 0.0 ----

# make the movie
./binary_stars.py --mass1 1.0 --mass2 1.0 -e 0.0 -a 10.0
mkmovie.py -o binary_mratio=1_e=0.0 binary_star_mratio=1.00_e=0.00_????.png

# store one of the frames as the preview frame for the webpage and clean-up
cp binary_star_mratio=1.00_e=0.00_0150.png binary_mratio=1_e=0.0_preview.png
rm -f binary_star_mratio=1.00_e=0.00_????.png


# ---- M1/M2 = 2.0; e = 0.0 ----

# make the movie
./binary_stars.py --mass1 2.0 --mass2 1.0 -e 0.0 -a 10.0
mkmovie.py -o binary_mratio=2_e=0.0 binary_star_mratio=2.00_e=0.00_????.png

# store one of the frames as the preview frame for the webpage and clean-up
cp binary_star_mratio=2.00_e=0.00_0150.png binary_mratio=2_e=0.0_preview.png
rm -f binary_star_mratio=2.00_e=0.00_????.png


# ---- M1/M2 = 1.0; e = 0.4 ----

# make the movie
./binary_stars.py --mass1 1.0 --mass2 1.0 -e 0.4 -a 10.0
mkmovie.py -o binary_mratio=1_e=0.4 binary_star_mratio=1.00_e=0.40_????.png

# store one of the frames as the preview frame for the webpage and clean-up
cp binary_star_mratio=1.00_e=0.40_0150.png binary_mratio=1_e=0.4_preview.png
rm -f binary_star_mratio=1.00_e=0.40_????.png


# ---- M1/M2 = 2.0; e = 0.4 ----

# make the movie
./binary_stars.py --mass1 2.0 --mass2 1.0 -e 0.4 -a 10.0
mkmovie.py -o binary_mratio=2_e=0.4 binary_star_mratio=2.00_e=0.40_????.png

# store one of the frames as the preview frame for the webpage and clean-up
cp binary_star_mratio=2.00_e=0.40_0150.png binary_mratio=2_e=0.4_preview.png
rm -f binary_star_mratio=2.00_e=0.40_????.png


# ---- M1/M2 = 4.0; e = 0.4 w/ energy ----

# make the movie
./binary_stars.py --mass1 4.0 --mass2 1.0 -e 0.4 -a 10.0 --annotate
mkmovie.py -o binary_mratio=4_e=0.4_energy binary_star_mratio=4.00_e=0.40_????_energy.png

# store one of the frames as the preview frame for the webpage and clean-up
cp binary_star_mratio=4.00_e=0.40_0150_energy.png binary_mratio=4_e=0.4_energy_preview.png
rm -f binary_star_mratio=4.00_e=0.40_????_energy.png

