#!/usr/bin/env python3
# Test function
# Used by the bin/example_*_.ini
# Yang, 2023.11.01
from numpy import loadtxt, sin, cos, random
from time import sleep

# read input from inp.dat
data = loadtxt(r"TestFunction_input.dat")

f = sin(data[0])**2 + cos(data[1])**2
sleep(0.5+0.1*random.random())

#open("TestFunction_output.dat",'w').write(str(f))
open("TestFunction_output.dat",'w').write("test "+str(f))
#open("TestFunction_output.dat",'w').write(
#"201  201# %s\n"%str(f)
#+ "BLOCK abc\n"
#+ "  0 0 %s"%str(f)
#)
