#!/usr/bin/env python

import os
import argparse
import re

# routines for natural/human sorting
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    # list.sort(key=natural_keys) will sort according 
    # to human sorting (or natural sorting)
    return [ atoi(c) for c in re.split('(\d+)', text) ]


parser = argparse.ArgumentParser(description="python driver to mencoder.  Take a list of files and make avi and mp4 movies.",
                                 usage="./mkmovie.py [options] <list of PNG files>",
                                 epilog="Note: for the mp4 movies, we duplicate the frames at the start to get\n"
                                        "around a bug (?) in mencoder that skips some frames at the start.\n\n"
                                        "M. Zingale (2013-02-19)\n\nR. Orvedahl (2014-07-18)\n\t-add more options\n\t"
                                        "-change calling sequences\n\n C. Gobat (2023-07-19)\n\t-use argparse instead of getopt",
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-f", "--fps", metavar="x", default=15, type=int, help="set frames per second to be x (default: 15)")
duplication = parser.add_mutually_exclusive_group()
duplication.add_argument("--double", action="store_true", help="double each frame (effectively 2x slower)")
duplication.add_argument("-N", dest="ncopies", metavar="num", default=1, type=int, help="create N copies of each frame (to really slow it down)")
parser.add_argument("-o", dest="prefix", metavar="prefix", default="movie", help="prefix for the movies (default: movie)")
movie_type = parser.add_mutually_exclusive_group()
movie_type.add_argument("--avi", dest="mp4", action="store_false", help="only generate an avi movie")
movie_type.add_argument("--mp4", dest="avi", action="store_false", help="only generate an mp4 movie")
parser.add_argument("--endframes", dest="endFrames", metavar="N", default=1, type=int, help="add N duplicate frames to the end, to make the last frame persist")
parser.add_argument("input_frames", nargs="+", help="list of image file paths/names to use as frames in the movie")

args = parser.parse_args()
frames = list(args.input_frames)

print("number of frames: ", len(frames))

# sort frames
frames.sort(key=natural_keys)

# create temporary files that list the movie frames -- these are inputs
# to mencoder

if (args.avi):
    f = open("_mkmovie1.list", "w")

    for img in frames:
        f.write("%s\n" % (img))
        if args.double:
            f.write("%s\n" % (img))
        elif args.ncopies > 1:
            for i in range(args.ncopies):
                f.write("%s\n" % (img))

    if args.endFrames > 1:
        n = 0
        while (n < args.endFrames-1):
            f.write("%s\n" % (frames[len(frames)-1]))        
            n += 1

    f.close()


if (args.mp4):
    # for mp4, we want some extra frames at the start
    f = open("_mkmovie2.list", "w")

    n = 0
    while (n < 28):
        f.write("%s\n" % (frames[0]))
        n += 1

    for img in frames:
        f.write("%s\n" % (img))
        if args.double:
            f.write("%s\n" % (img))
        elif args.ncopies > 1:
            for i in range(args.ncopies):
                f.write("%s\n" % (img))

    if args.endFrames > 1:
        n = 0
        while (n < args.endFrames-1):
            f.write("%s\n" % (frames[len(frames)-1]))        
            n += 1

    f.close()



# make the movies
if (args.avi):
    str1 = "mencoder mf://@_mkmovie1.list -ovc lavc -lavcopts "
    str2 = "vcodec=msmpeg4v2:vbitrate=3000:vhq:mbd=2:trell "
    str3 = "-mf type=png:fps="+str(args.fps)
    str4 = f" -o {args.prefix}.avi"
    os.system(str1+str2+str3+str4)

if (args.mp4):
    str1 = "mencoder mf://@_mkmovie2.list -of lavf -lavfopts format=mp4 "
    str2 = "-ss 1 -ovc x264 -x264encopts crf=20.0:nocabac:level_idc=30:"
    str3 = "global_header:threads=2 -fps "+str(args.fps)
    str4 = f" -o {args.prefix}.mp4"
    os.system(str1+str2+str3+str4)

    str1 = "mencoder mf://@_mkmovie2.list -ovc x264 -x264encopts crf=10:"
    str2 = "me=umh:subq=9:nr=100:global_header -of lavf -lavfopts format=mp4 "
    str3 = "-fps "+str(args.fps)
    str4 = f" -o {args.prefix}_hg.mp4"
    os.system(str1+str2+str3+str4)


# remove the files
if (args.avi):
    os.remove("_mkmovie1.list")
if (args.mp4):
    os.remove("_mkmovie2.list")

if (args.avi):
    print(f"\n Avi Movie: {args.prefix}.avi")
if (args.mp4):
    print(f"\n Mp4 Movie: {args.prefix}.mp4")
    print(f"\n Mp4 Movie: {args.prefix}_hg.mp4")

