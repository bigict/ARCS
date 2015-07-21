#!/usr/bin/python

Usage = '''
Author : Zheng quangang

$python links_dis.py -D <workspace>
                     -L <links_dis.data>
                     -U <UR>
'''
import sys
import getopt

opts, args = getopt.getopt(sys.argv[1:],"D:L:U:")
if len(opts) != 3:
    print Usage
    sys.exit()

for op,value in opts:
    if op == "-D":
        prefix = value
    elif op == "-L":
        ld = value
    else:
        ur = value

links = {}
num = {}
lens = []
flag = []

cin1 = open(ld,"r")
cin2 = open(ur,"r")
try:
    while 1:
        line = cin1.next().strip()
        array = line.split()
        i = int(array[0])
        j = int(array[1])
        if (i,j) in links:
            links[(i,j)].append(int(array[2]))
            num[(i,j)] += 1
        else:
            links[(i,j)] = []
            links[(i,j)].append(int(array[2]))
            num[(i,j)] = 1
except:
    info = sys.exc_info()
    print info[0],info[1]

try:
    while 1:
        line = cin2.next().strip()
        array = line.split()
        i = int(array[0])
        lens.append(int(array[1]))
        flag.append(array[2])
except:
    info = sys.exc_info()
    print info[0],info[1]

cin1.close()
cin2.close()


cout = open(prefix+"/Cyto_links_dis2","w")
cout.write("i\tlen(i)\tj\tlen(j)\tUU/UR\tdis...\n")
for (i,j) in links:
    cout.write(str(len(links[(i,j)]))+"\t"+str(i)+"\t"+str(lens[i])+"\t"+str(j)+"\t"+str(lens[j])+"\t"+flag[i]+flag[j]+"\t")
    for k in links[(i,j)]:
        cout.write(str(k)+",")
    cout.write("\n")
cout.close()







