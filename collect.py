#!/bin/python3

import os
import sys

ref = f'''
    Usage:
        {__file__} <out dir> <file name> <result file>

    <out dir>:
        The path to the directory, contains the output of MPI code.

    <file name>:
        The name of the files that should be collected among the
        subdirectories.

    <result file>:
        The name of the file, that will be used for storing the
        result of the collecting.
'''

if (sys.argv[1] == "-h" or sys.argv[1] == "--help"):
    print(ref)
    sys.exit()

directory = sys.argv[1]
filename = sys.argv[2]
output_file = sys.argv[3]

coords = [
    (int(sub.split('-')[0]), int(sub.split('-')[1]))
    for sub in os.listdir(directory)
]

x_count = max(coords, key = lambda item: item[0])[0] + 1
y_count = max(coords, key = lambda item: item[1])[1] + 1

with open(output_file, "w") as out:
    for y in range(y_count):
        files = [
            open(os.path.join(directory, f"{x}-{y}", filename), "r")
            for x in range(x_count)
        ]

        for strs in zip(*files):
            out.write(' '.join(map(str.strip, strs)) + '\n')

        for f in files:
            f.close()