#!/usr/bin/python
#
#import sys
#from os import listdir
#from os.path import isfile, join
#import FortranFile as ff
#import re
#
#if(len(sys.argv) != 3):
#    print("ERROR wrong number of arguments.")
#    sys.exit(1)
#
## Keep the mseek'th timestep and truncate everything else
#mseek = int(sys.argv[2])
#mypath = sys.argv[1]
#
#files = [ff.FortranFile(join(mypath,f), mode='r+') for f in listdir(mypath) 
#                    if isfile(join(mypath,f)) 
#                    and 'outflowing' not in f]
#
#prev_m = -1
#for f in files:
#    index, trunc = f.index(repair=True)
#    if trunc:
#        print('{} was repaired'.format(f.name))
#    # Timestep numbers are the even entries in index.
#    # Iterate over those, but not the last one since that
#    # is not the start of a record, it is the end of the file.
#    for i in index[::2][:-1]:
#        f.seek(i)
#        m, = f.readInts()
#        if m < 0:
#            raise ValueError("Timestep numbers should always be positive.")
#        if m < prev_m:
#            raise ValueError("Timestep numbers should be strictly increasing")
#
#        if m > mseek:
#            # go back to before this step
#            f.seek(i)
#            # Cut off everything that is going to be redone anyway
#            f.truncate()
#            # Then we're done with this file
#            break
#
#        m_prev = m
