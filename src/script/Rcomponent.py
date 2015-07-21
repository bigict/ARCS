#!/usr/bin/python

import sys
import os

if len(sys.argv) != 3:
	print "Usage: Rcomponent.py <ITER> <workspace> <threshold of len>"
	os._exit(0)
threshold = 50
if len(sys.argv) == 4:
    threshold = int(sys.argv[3])

print "[ Checkpoint ] Building new component using repeat contigs...\n";
cdbgl = []
R = []
index = 0
fin = open(sys.argv[2]+"/cdbg_copy_number.fa", 'r')
try:
    while 1:
        line = fin.next().strip()
        array = line.split()
        line = fin.next().strip()
        if int(array[1]) > 1 and int(array[1]) <= 10:
            R.append(index)
        cdbgl.append(len(line))
        index += 1
except:
    info = sys.exc_info()
    print info[0],info[1]
fin.close()

count = 0
fout = open(sys.argv[2]+"/component_"+str(10),"w")
fin = open(sys.argv[2]+"/component_"+str(1), 'r')
try:
    while 1:
        line = fin.next().strip()
        line = fin.next().strip()
        array = line.split()
        if len(array) == 1:
            if cdbgl[int(array[0])] < threshold:
                line = fin.next()
                continue
        fout.write(">component "+str(count)+"\n")
        fout.write(line+"\n")
        line = fin.next().strip()
        fout.write(line+"\n")
        count += 1
except:
    info = sys.exc_info()
    print info[0],info[1]
fin.close()

for i in R:
    if cdbgl[i] >= threshold:
        fout.write(">component "+str(count)+"\n")
        fout.write(str(i)+"\n")
        fout.write("\n")
        count += 1

fout.close()
print "[ Done ] "
