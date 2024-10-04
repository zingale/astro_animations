ffmpeg -framerate 15 -f image2 -pattern_type glob -i "summer_solstice*.png" -vcodec mpeg4 -c:v libx264 -crf 20 summer_solstice.mp4
ffmpeg -framerate 15 -f image2 -pattern_type glob -i "winter_solstice*.png" -vcodec mpeg4 -c:v libx264 -crf 20 winter_solstice.mp4
ffmpeg -framerate 15 -f image2 -pattern_type glob -i "equinox*.png" -vcodec mpeg4 -c:v libx264 -crf 20 equinox.mp4
