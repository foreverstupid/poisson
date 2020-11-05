#!/bin/python3

import matplotlib as mpl
import matplotlib.pyplot as plt
import csv
import sys
import numpy as np

x1 = 0.0
x2 = 4.0
y1 = 0.0
y2 = 3.0

xcount = int(sys.argv[1])
ycount = int(sys.argv[2])
out_file = sys.argv[3]
plot_file_name = sys.argv[4]

xrange = np.linspace(x1, x2, xcount)
yrange = np.linspace(y1, y2, ycount)
xs, ys = np.meshgrid(xrange, yrange)

with open(out_file, "r") as f:
    reader = csv.reader(f, delimiter = ' ')
    zs = np.array([
        [ float(val) for val in item ] for item in reader
    ])

fig = plt.figure()
hm_axis = fig.add_subplot(1, 1, 1, projection="3d")

real_zs = np.sqrt(xs * ys + 4)
hm = hm_axis.plot_surface(xs, ys, zs)
fig.savefig(plot_file_name)