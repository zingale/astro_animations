ffmpeg -framerate 15 -f image2 -pattern_type glob -i "esc*.png" -vcodec mpeg4 -c:v libx264 -crf 20 escape.mp4
