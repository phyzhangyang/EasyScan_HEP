#!/usr/bin/env python3
# Ackley function
# Yang, 2019.01.23
from numpy import loadtxt, sqrt, exp, cos, pi

# read input from inp.dat
data = loadtxt(r"AckleyFunction_input.dat")

x = data[0]
y = data[1]

f = -20*exp( -0.2*sqrt(0.5*(x*x + y*y)) ) - exp( 0.5*(cos(2*pi*x) + cos(2*pi*y)) ) + 22.718282;


open("AckleyFunction_output.dat",'w').write(str(f))
