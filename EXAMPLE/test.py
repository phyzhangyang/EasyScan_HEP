#!/usr/bin/env python
# For test ES
# by ZY
from numpy import *

data = loadtxt(r"input.txt")

x = data[0]
y = data[1]

z = sin(x)**2 + cos(y)**2

open("output.txt",'w').write(str(z))