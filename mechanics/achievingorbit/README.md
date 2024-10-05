ffmpeg -framerate 30 -f image2 -pattern_type glob -i "ach*.png" -vcodec mpeg4 -c:v libx264 -crf 20 -pix_fmt yuv420p -movflags +faststart achievingorbit.mp4
