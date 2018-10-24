#! /bin/python3

import argparse
import subprocess
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

subprocess.run(["mkdir", "-p", "heatmap_frames"], capture_output=True, text=True);

# yeah, this way we have to keep it all in memory, but if it becomes a problem I'll figure something
# else out
print("Calculating transient temperature distribution...")
calculations = subprocess.run(["./" + args.exec, str(args.npts), str(args.nt), str(args.dt),
    str(args.alpha)], capture_output=True, text=True)

frames = list(filter(None, calculations.stdout.split('\n')))
# format: [node 0,0],[node 1,0],...;[node 0,1],[node 1,1],...;...,[node (num_points - 1),(num_points - 1)];\n
frame_num = 0
print("Generating animation frames...")
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
    plot.imsave("heatmap_frames/" + str(frame_num), frame_np, cmap='hot');
    if frame_num % 100 == 0:
        print("\t" + str(frame_num) + " frames generated...")
    frame_num = frame_num + 1

print("Starting animation conversion...");
cmd = "ffmpeg -y -i heatmap_frames/%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p -vf scale=-1:480 -framerate " + str(args.framerate) + " " + args.outfile
conversion = subprocess.run(cmd.split(' '), capture_output=True, text=True)
frame_del = subprocess.run(["rm", "-rf", "heatmap_frames"], capture_output=True, text=True)
print("...done! " + args.outfile + " generated")
