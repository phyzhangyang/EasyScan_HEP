#!/usr/bin/env python
# Test function
# Used by the bin/example.ini
# Yang, 2019.01.23
from numpy import loadtxt, sin, cos

# read input from inp.dat
data = loadtxt(r"TestFunction_input.dat")

f = sin(data[0])**2 + cos(data[1])**2

open("TestFunction_output.dat",'w').write("test "+str(f))
#open("TestFunction_output.dat",'w').write(
#"201  201# %s\n"%str(f)
#+ "BLOCK abc\n"
#+ "  0 0 %s"%str(f)
#)
