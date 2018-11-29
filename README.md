# Installation
1. clone (or otherwise download) this repo
2. `make`
   * this generates the actual calculator, written in C
   * the only dependency should be `gcc`
   * you can run the calculator on its own, but you're just going to see a mess of numbers on stdout

# Usage
1. `python heatmap.py [output_file.mp4]`
   * this wrapper invokes the calculator and generates an animation
   * make sure you have python3, ffmpeg, numpy, and matplotlib installed
      * use your distro's package installer for the first 2 (e.g. `sudo pacman -S python3 ffmpeg`)
      * try `sudo pip install numpy matplotlib` for the last 2
   * see advanced options with `./heatmap.py -h`
   * you don't have to use mp4; any container that supports h264 should work (ffmpeg is magical)
2. profit!

# Notes
* The threaded implementation actually works slower on my machine - the time saved by parallel
  calculation of the interior is smaller than the overhead of spinning up threads. As the number of
  nodal points increases and each time step takes longer to calculate, the threads will show a
  performance gain, but the total runtime will probably be unusably high by the time the threads
  save any time.
