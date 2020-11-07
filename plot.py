#!/bin/python3

import matplotlib as mpl
import matplotlib.pyplot as plt
import csv
import sys
import numpy as np

ref = f'''
    Usage:
        {__file__} <data file> <out file> [-d]

    <data file>:
        The file, contains the data.

    <out file>:
        The file, that will be used as an output one.

    -d:
        If specified then it plots the difference between the
        given and expected data.
'''

if (sys.argv[1] == "-h" or sys.argv[1] == "--help"):
    print(ref)
    sys.exit()

x1 = 0.0
x2 = 4.0
y1 = 0.0
y2 = 3.0

data_file = sys.argv[1]
plot_file_name = sys.argv[2]

with open(data_file, "r") as f:
    reader = csv.reader(f, delimiter = ' ')
    zs = np.array([
        [ float(val) for val in item ] for item in reader
    ])

xrange = np.linspace(x1, x2, zs.shape[1])
yrange = np.linspace(y1, y2, zs.shape[0])
xs, ys = np.meshgrid(xrange, yrange)

fig = plt.figure()
hm_axis = fig.add_subplot(1, 1, 1, projection="3d")

if (len(sys.argv) > 3 and sys.argv[3] == "-d"):
    zs = zs - np.sqrt(xs * ys + 4)

hm = hm_axis.plot_surface(xs, ys, zs)
fig.savefig(plot_file_name)