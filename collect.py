#!/bin/python3

import os
import sys

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