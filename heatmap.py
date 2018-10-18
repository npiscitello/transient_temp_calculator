#! /bin/python3

import argparse
import subprocess
import os
import numpy as np
import matplotlib.pyplot as plot

parser = argparse.ArgumentParser()
parser.add_argument("outfile", help="filename of the generated animation")
parser.add_argument("-r", "--framerate", type=int, help="frames per second of the generated animation", default=100)
parser.add_argument("-p", "--npts", type=int, help="number of nodes in x and y", default=40)
parser.add_argument("-n", "--nt", type=int, help="number of time steps to take", default=1000)
parser.add_argument("-d", "--dt", type=float, help="size of the time steps, in seconds", default=0.1)
parser.add_argument("-a", "--alpha", type=float, help="thermal diffusivity, in m^2/s", default=0.001)
parser.add_argument("-e", "--exec", help="executable name of the C calculator", default="transient")
args = parser.parse_args()

os.makedirs("heatmap_frames")

# yeah, this way we have to keep it all in memory, but if it becomes a problem I'll figure something
# else out
calculations = subprocess.run(["./" + args.exec, str(args.npts), str(args.nt), str(args.dt),
    str(args.alpha)], capture_output=True, text=True)

frames = list(filter(None, calculations.stdout.split('\n')))
# format: [node 0,0],[node 1,0],...;[node 0,1],[node 1,1],...;...,[node (num_points - 1),(num_points - 1)];\n
frame_num = 0
for frame_raw in frames:
    rows = list(filter(None, frame_raw.split(';')))
    frame = []
    for row_raw in rows:
        vals = list(filter(None, row_raw.split(',')))
        row = []
        for val in vals:
            row.append(float(val))
        frame.append(row)
    frame_np = np.array(frame)
    plot.imsave("heatmap_frames/" + args.outfile + "_" + str(frame_num), frame_np, cmap='hot');
    frame_num = frame_num + 1

#ffmpeg -framerate 100 -i heatmap_frames/test_%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p test.mp4
