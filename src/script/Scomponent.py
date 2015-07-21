#!/usr/bin/python

import sys
import os

if len(sys.argv) < 2:
	print "Usage: Scomponent.py <ITER> <workspace> <scaf_seq_with_gaps> <threshold of len> "
	os._exit(0)
threshold = 50
if len(sys.argv) == 5:
    threshold = int(sys.argv[4])

print "[Checkpoint]\tBuilding new component and cdbg based on mer.sca_seq_with_gaps\n";
cdbgl = []
index = 0
fin = open(sys.argv[2]+"/"+sys.argv[3], 'r')
count = 0
fout = open(sys.argv[2]+"/component_"+str(20),"w")
cout = open(sys.argv[2]+"/cdbg_sca.fa","w")
try:
    while 1:
        line = fin.next().strip()
        line = fin.next().strip()
        if len(line) < threshold:
            continue
        fout.write(">component_"+str(len(line))+" "+str(count)+"\n")
        fout.write(str(count)+"\n\n")

        cout.write(">seq_"+str(count)+"\t"+str(1)+"\n")
        cout.write(line+"\n")
        count += 1
except:
    info = sys.exc_info()
    print info[0],info[1]
fin.close()
fout.close()
cout.close()

print "[ Done ] "
