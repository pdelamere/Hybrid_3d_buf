#!/usr/bin/python3

import sys

if(len(sys.argv) != 3):
    print("ERROR wrong number of arguments.")
    sys.exit(1)

mseek = int(sys.argv[2])
mypath = sys.argv[1]

from os import listdir
from os.path import isfile, join

files = [open(join(mypath,f),'r+b') for f in listdir(mypath) if isfile(join(mypath,f))]

print(list(map(lambda x:x.name,files)))
print('waiting...')
input()

import struct

for f in files:
    r = f.read(4)
    if (r==b''):
        continue

    f.seek(0)
    m=0
    while(m < mseek):
        ## Read first record size
        r = f.read(4)
        assert (r != b''), "m="+str(m)+" mseek="+str(mseek)+" file="+f.name
        rl1, = struct.unpack('i',r)
        assert(rl1==4), "m="+str(m)+" mseek="+str(mseek)+" file="+f.name

        ## Read step number
        r = f.read(4)
        assert(r != b''), "m="+str(m)+" mseek="+str(mseek)+" file="+f.name
        m, = struct.unpack('i',r)

        ## Read second record size
        r = f.read(4)
        assert(r != b''), "m="+str(m)+" mseek="+str(mseek)+" file="+f.name
        rl2, = struct.unpack('i',r)
        assert(rl2==4), "m="+str(m)+" mseek="+str(mseek)+" file="+f.name
        
        ## Read first record size for the data record
        r = f.read(4)
        assert(r != b''), "m="+str(m)+" mseek="+str(mseek)+" file="+f.name
        rl1, = struct.unpack('i',r)

        ## Skip over the data
        f.seek(rl1,1)

        ## Read second record size
        r = f.read(4)
        assert(r != b''), "m="+str(m)+" mseek="+str(mseek)+" file="+f.name
        rl2, = struct.unpack('i',r)
        assert(rl1==rl2), "m="+str(m)+" mseek="+str(mseek)+" file="+f.name
    ## cut off everything that's going to be redone anyway
    f.truncate()
    f.close()
