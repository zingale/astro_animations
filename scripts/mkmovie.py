#!/usr/bin/env python

desc = """ 
python driver to mencoder.  Take a list of files and make avi and mp4 movies.

options:

  --double: double each frame (effectively 2x slower)

  -N num: create N copies of each frame (to really slow it down)
 
  -o prefix: prefix for the movies (default: movie)

  --endframes N: add N duplicate frames to the end, to make the last 
                 frame persist

Note: for the mp4 movies, we duplicate the frames at the start to get
around a bug (?) in mencoder that skips some frames at the start.

M. Zingale (2013-02-19)

"""

import sys
import os
import getopt


if len(sys.argv) == 1:
    print desc
    sys.exit(2)


prefix = "movie"
endFrames = 1
double = 0
ncopies = 1

try: opts, next = getopt.getopt(sys.argv[1:], "o:N:", ["double", "endframes="])
except getopt.GetoptError:
    sys.exit("invalid calling sequence")

for o, a in opts:

    if o == "-o":
        prefix = a

    elif o == "--double":
        double = 1

    elif o == "-N":
        ncopies = int(a)

    elif o == "--endframes":
        endFrames = int(a)

try: frames = next[0:]
except IndexError:
    sys.exit("no frames specified")

print "number of frames: ", len(frames)


# create temporary files that list the movie frames -- these are inputs
# to mencoder

f = open("_mkmovie1.list", "w")

for img in frames:
    f.write("%s\n" % (img))
    if double:
        f.write("%s\n" % (img))
    elif ncopies > 1:
        for i in range(ncopies):
            f.write("%s\n" % (img))

if endFrames > 1:
    n = 0
    while (n < endFrames-1):
        f.write("%s\n" % (frames[len(frames)-1]))        
        n += 1

f.close()


# for mp4, we want some extra frames at the start
f = open("_mkmovie2.list", "w")

n = 0
while (n < 28):
    f.write("%s\n" % (frames[0]))
    n += 1

for img in frames:
    f.write("%s\n" % (img))
    if double:
        f.write("%s\n" % (img))
    elif ncopies > 1:
        for i in range(ncopies):
            f.write("%s\n" % (img))

if endFrames > 1:
    n = 0
    while (n < endFrames-1):
        f.write("%s\n" % (frames[len(frames)-1]))        
        n += 1

f.close()



# make the movies
os.system("mencoder mf://@%s -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=3000:vhq:mbd=2:trell -mf type=png -o %s.avi" % ("_mkmovie1.list", prefix) )

os.system("mencoder mf://@%s -of lavf -lavfopts format=mp4 -ss 1 -ovc x264 -x264encopts crf=20.0:nocabac:level_idc=30:global_header:threads=2 -o %s.mp4" % ("_mkmovie2.list", prefix) )

os.system("mencoder mf://@%s -ovc x264 -x264encopts crf=10:me=umh:subq=9:nr=100:global_header -of lavf -lavfopts format=mp4 -o %s_hg.mp4" % ("_mkmovie2.list", prefix) )



# remove the files
os.remove("_mkmovie1.list")
os.remove("_mkmovie2.list")




