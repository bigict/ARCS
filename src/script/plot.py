#!/usr/bin/python

import sys
import os
import getopt


fil = ""
opts, var = getopt.getopt(sys.argv[1:], "i:D:F:")
if len(opts) != 2:
	print "Usage: plot.py -D <workspace> -F <file>"
	os._exit(0)
for op,value in opts:
	if op == '-i':
		ite = int (value)
	elif op == '-D':
		prefix = value
	elif op == '-F':
		fil = value

index = 0
fin = open(fil, 'r')
table = {}
try:
    while 1:
        line1 = fin.next().strip()
        name = line1[5:8]
        line2 = fin.next().strip()
        array = line2.split()
        #print array
        index = -100
        fout0 = open("tmp", 'w')
        t = []
        for i in array:
            fout0.write(str(index) + "\t" + i + "\n")
            t.append(i)
            index += 1
        fout0.close()
        table[name] = t
        cmd = "./plot.sh " + prefix+name+".pdf" + " tmp " + name 
        print cmd
        print os.system(cmd)
except:
    info = sys.exc_info()
    print info[0],info[1]
fin.close()
'''
fout = open("tmp","w")
key = []
for i in table:
    if int(i[0]) != 0 or int(i[2]) != 0:
        key.append(i)
for i in key:
    del table[i]
index = -100
lens = -1
for i in table:
    if lens > 0 and len(table[i]) >= lens:
        continue
    else:
        lens = len(table[i])
key = []
for i in table:
    print i
    key.append(i)
for i in range(lens):
    fout.write(str(index)+"\t")
    for i in key:
        fout.write(i+"\t"+str(table[i][index+100])+"\t")
    index += 1
    fout.write("\n")
fout0.close()
cmd = "./plot3.sh " + prefix+"0x0.pdf" + " tmp " + name 
print cmd
print os.system(cmd)
'''

print "[Done] "
