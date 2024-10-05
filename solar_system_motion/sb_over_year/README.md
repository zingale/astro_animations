ffmpeg -framerate 15 -f image2 -pattern_type glob -i "stonybrook_noon_moll_*.png" -vcodec mpeg4 -c:v libx264 -crf 20 -pix_fmt yuv420p -movflags +faststart stonybrook_noon_year_moll.mp4

ffmpeg -framerate 15 -f image2 -pattern_type glob -i "stonybrook_noon_ortho_*.png" -vcodec mpeg4 -c:v libx264 -crf 20 -pix_fmt yuv420p -movflags +faststart stonybrook_noon_year.mp4


