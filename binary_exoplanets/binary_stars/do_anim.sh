#!/bin/bash

# ---- M1/M2 = 1.0; e = 0.0 ----

# make the movie
./binary_stars.py --mass1 1.0 --mass2 1.0 -e 0.0 -a 10.0
ffmpeg -framerate 15 -f image2 -pattern_type glob -i "binary_star_mratio=1.00_e=0.00_????.png" -vcodec mpeg4 -c:v libx264 -crf 20 -pix_fmt yuv420p -movflags +faststart binary_mratio=1_e=0.0.mp4

rm -f binary_star_mratio=1.00_e=0.00_????.png


# ---- M1/M2 = 2.0; e = 0.0 ----

# make the movie
./binary_stars.py --mass1 2.0 --mass2 1.0 -e 0.0 -a 10.0
ffmpeg -framerate 15 -f image2 -pattern_type glob -i "binary_star_mratio=2.00_e=0.00_????.png" -vcodec mpeg4 -c:v libx264 -crf 20 -pix_fmt yuv420p -movflags +faststart binary_mratio=2_e=0.0.mp4

rm -f binary_star_mratio=2.00_e=0.00_????.png


# ---- M1/M2 = 1.0; e = 0.4 ----

# make the movie
./binary_stars.py --mass1 1.0 --mass2 1.0 -e 0.4 -a 10.0
ffmpeg -framerate 15 -f image2 -pattern_type glob -i "binary_star_mratio=1.00_e=0.40_????.png" -vcodec mpeg4 -c:v libx264 -crf 20 -pix_fmt yuv420p -movflags +faststart binary_mratio=1_e=0.4.mp4

rm -f binary_star_mratio=1.00_e=0.40_????.png


# ---- M1/M2 = 2.0; e = 0.4 ----

# make the movie
./binary_stars.py --mass1 2.0 --mass2 1.0 -e 0.4 -a 10.0
ffmpeg -framerate 15 -f image2 -pattern_type glob -i "binary_star_mratio=2.00_e=0.40_????.png" -vcodec mpeg4 -c:v libx264 -crf 20 -pix_fmt yuv420p -movflags +faststart binary_mratio=2_e=0.4.mp4

rm -f binary_star_mratio=2.00_e=0.40_????.png


# ---- M1/M2 = 4.0; e = 0.4 w/ energy ----

# make the movie
./binary_stars.py --mass1 4.0 --mass2 1.0 -e 0.4 -a 10.0 --annotate
ffmpeg -framerate 15 -f image2 -pattern_type glob -i "binary_star_mratio=4.00_e=0.40_????_energy.png" -vcodec mpeg4 -c:v libx264 -crf 20 -pix_fmt yuv420p -movflags +faststart binary_mratio=4_e=0.4_energy.mp4

rm -f binary_star_mratio=4.00_e=0.40_????_energy.png

