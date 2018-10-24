# Installation
1. clone (or otherwise download) this repo
2. `make`
  * this generates the actual calculator, written in C
  * the only dependency should be `gcc`
  * you can run the calculator on its own, but you're just going to see a mess of numbers on stdout

# Usage
1. `./heatmap.py [output_file.mp4]`
  * this wrapper invokes the calculator and generates an animation
  * make sure you have python3, ffmpeg, numpy, and matplotlib installed
    * use your distro's package installer for the first 2 (e.g. `sudo pacman -S python3 ffmpeg`)
    * try `sudo pip install numpy matplotlib` for the last 2
  * see advanced options with `./heatmap.py -h`
  * you don't have to use mp4; any container that supports h264 should work (ffmpeg is magical)
2. profit!
