#!/usr/bin/env python

import sys
import numpy as np
import re

if len(sys.argv) != 4:
    print("Usage: %s rle_file startx starty.")
    sys.exit(0)

f = open(sys.argv[1], 'r')

x = 0
y = 0
r = ''
data = ''
for l in f:
    if l[0] == '#':
        continue
    if l[0] == 'x':
        v = l.split(',')
        x = int(v[0].strip().split('=')[1].strip())
        y = int(v[1].strip().split('=')[1].strip())
        r = v[2].strip().split('=')[1].strip()
        continue
    if l.strip('') != '':
        data = '%s%s' % (data, l.strip())

print(x, y, r)
print(data.split('$'))

a = np.zeros((x, y), dtype = int)

def prline(a):
    s = re.split('([bo])', a)
    l = ''
    mul = 1
    for i in s:
        if i == 'b' or i == 'o':
            l = '%s%s' % (l, i * mul)
        else:
            if i == '':
                mul = 1
            else:
                mul = int(i)
    return l
j = 0
for l in data.split('$'):
    d = prline(l.strip('!'))
    for i in range(len(d)):
        if d[i] == 'o':
            a[i,j] = 1
    j = j + 1

stx = int(sys.argv[2])
sty = int(sys.argv[3])
for i in range(x):
    for j in range(y):
        if a[i,j] == 1:
            print(i+stx,j+sty)