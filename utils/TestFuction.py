#!/usr/bin/env python
# Test function
# Used by the bin/example.ini
# Yang, 2019.01.23
from numpy import loadtxt, sin, cos

# read input from inp.dat
data = loadtxt(r"TestFuction_input.dat")

f = sin(data[0])**2 + cos(data[1])**2

open("TestFuction_output.dat",'w').write("test "+str(f))
