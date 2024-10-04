ffmpeg -framerate 30 -f image2 -pattern_type glob -i "chan*.png" -vcodec mpeg4 -c:v libx264 -crf 20 changing_orbit.mp4

