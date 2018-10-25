#!/usr/bin/env python
import linecache
import sys
import os

if len(sys.argv) == 4 :
    if os.path.exists(sys.argv[2])==True: os.remove(sys.argv[2]) 
    print int(float(open(sys.argv[3]).readline().rstrip()))
    try:
        open(sys.argv[2],'w').write(linecache.getline(sys.argv[1],int(float(open(sys.argv[3]).readline().rstrip()))))
    except:
        print "Error: 'read_one_line.py' need 3 input parameters: input_file_name, output_file_name, and line number."
else:
    print "Error: 'read_one_line.py' need 3 input parameters: input_file_name, output_file_name, and line number."


