ffmpeg -framerate 0.5 -f image2 -pattern_type glob -i "*.png" -vcodec mpeg4 -c:v libx264 -crf 20 -pix_fmt yuv420p -movflags +faststart decay.mp4
