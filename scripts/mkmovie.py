#!/usr/bin/env python

import os
import argparse
import re
import subprocess

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
                                        "M. Zingale (2013-02-19)\n\nR. Orvedahl (2014-09-02)\n\t-add more options\n\t"
                                        "-change calling sequences\n\n C. Gobat (2023-07-19)\n\t-use argparse instead of getopt;"
                                        " subprocess instead of os.system\n\t-add OpenCV fallback if mencoder fails\n\t"
                                        "-general modernization and cleanup",
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


if (args.avi):
    avi_frames = []

    for img in frames:
        avi_frames.append(img)
        if args.double:
            avi_frames.append(img)
        elif args.ncopies > 1:
            for i in range(args.ncopies):
                avi_frames.append(img)

    avi_frames += [frames[-1]]*(args.endFrames-1)

    with open("_mkmovie1.list", "w+") as f:
        # input list file for mencoder
        f.writelines(avi_frames)

    try:
        p = subprocess.run("mencoder mf://@_mkmovie1.list -ovc lavc -lavcopts "
                           "vcodec=msmpeg4v2:vbitrate=3000:vhq:mbd=2:trell "
                           f"-mf type=png:fps={args.fps}"
                           f" -o {args.prefix}.avi".split(),
                           check=True, )
    except:
        import cv2
        first_img = cv2.imread(frames[0])
        height, width, layers = first_img.shape # determine frame dimensions
        avi_writer = cv2.VideoWriter(args.prefix+".avi", cv2.VideoWriter.fourcc(*"MJPG"),
                                     args.fps, (width, height))
        for frame in avi_frames:
            avi_writer.write(cv2.imread(frame))
        cv2.destroyAllWindows()
        avi_writer.release()
    finally:
        os.remove("_mkmovie1.list")

if (args.mp4):
    
    # for mp4, we want some extra frames at the start
    mp4_frames = [frames[0]]*28

    for img in frames:
        mp4_frames.append(img)
        if args.double:
            mp4_frames.append(img)
        elif args.ncopies > 1:
            for i in range(args.ncopies):
                mp4_frames.append(img)

    mp4_frames += [frames[-1]]*(args.endFrames-1)
    
    with open("_mkmovie2.list", "w+") as f:
        # input list file for mencoder
        f.writelines(mp4_frames)

    try:
        subprocess.run("mencoder mf://@_mkmovie2.list -of lavf -lavfopts format=mp4 "
                       "-ss 1 -ovc x264 -x264encopts crf=20.0:nocabac:level_idc=30:"
                       f"global_header:threads=2 -fps {args.fps}"
                       f" -o {args.prefix}.mp4".split(), check=True)
    except:
        import cv2
        mp4_writer = cv2.VideoWriter(args.prefix+".mp4", cv2.VideoWriter.fourcc(*"mp4v"),
                                     args.fps, (width, height))
        for frame in mp4_frames:
            mp4_writer.write(cv2.imread(frame))
        cv2.destroyAllWindows()
        mp4_writer.release()
    try:
        subprocess.run("mencoder mf://@_mkmovie2.list -ovc x264 -x264encopts crf=10:"
                       "me=umh:subq=9:nr=100:global_header -of lavf -lavfopts format=mp4 "
                       f"-fps {args.fps} -o {args.prefix}_hg.mp4".split(), check=True)
    except:
        pass
    finally:
        os.remove("_mkmovie2.list")


if os.path.exists(args.prefix+".avi"):
    print(f"\n Avi Movie: {args.prefix}.avi")
if os.path.exists(args.prefix+".mp4"):
    print(f"\n Mp4 Movie: {args.prefix}.mp4")
if os.path.exists(args.prefix+"_hg.mp4"):
    print(f"\n Mp4 Movie: {args.prefix}_hg.mp4")
