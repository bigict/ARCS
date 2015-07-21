#! /usr/bin/python

Usage = '''
Author : Zheng QuanGang
Function : 

$python delCyto.py -D <workspace>
		   
'''

import sys,getopt
import os,math
import math

opts, args = getopt.getopt(sys.argv[1:],"I:K:C:S:D:M:R:")
if len(opts) < 0:
	print Usage
	sys.exit()
	
for op, value in opts:
	if op == '-K':
		K = int(value)
	elif op == '-C':
		cdbg = value
	elif op == '-S':
		scf = value
	elif op == '-D':
		prefix = value
	elif op == '-M':
		com = value
	elif op == '-R':
		ref = value
	elif op == '-I':
		I = int(value)
	else:
		print Usage
		sys.exit()
flag = {}
ctg_pos = {}
emdis = {}
refdis = {}
links = {}
lens = {}

cin = open(prefix+"/Cyto_contig_pos","r")
try:
    while 1:
        line = cin.next().strip()
        array = line.split()
        flag[array[0]] = array[1]
        ctg_pos[array[0]] = []
        listp = array[2].split(',')
        for i in listp:
            if i != '':
                ctg_pos[array[0]].append(int(i))
except:
    info = sys.exc_info()
    print info[0],info[1]
    cin.close()

cin = open(prefix+"/Cyto_emdis","r")
try:
    while 1:
        line = cin.next().strip()
        array = line.split()
        i = int(array[0])
        j = int(array[2])
        emdis[(i,j)] = []
        dis = array[4].split(',')
        if float(dis[0]) != 0.0:
            emdis[(i,j)].append(dis[0])
        if float(dis[1]) != 0.0:
            emdis[(i,j)].append(dis[1])
except:
    info = sys.exc_info()
    print info[0],info[1]
    cin.close()
	
cin = open(prefix+"/Cyto_refdis","r")
try:
    while 1:
        line = cin.next().strip()
        array = line.split()
        i = int(array[0])
        j = int(array[2])
        refdis[(i,j)] = []
        dis = array[4].split(',')
        if float(dis[0]) != 0.0:
            refdis[(i,j)].append(float(dis[0]))
        if float(dis[1]) != 0.0:
            refdis[(i,j)].append(float(dis[1]))
except:
    info = sys.exc_info()
    print info[0],info[1]
    cin.close()

cin = open(prefix+"/Cyto_links_dis","r")
try:
    while 1:
        line = cin.next().strip()
        array = line.split()
        i = int(array[0])
        j = int(array[1])
        links[(i,j)] = []
        dis = array[2].split(',')
        for d in dis:
            if d != '':
                links[(i,j)].append(int(d))
except:
    info = sys.exc_info()
    print info[0],info[1]
    cin.close()

cin = open(prefix+"/Cyto_node","r")
try:
    while 1:
        line = cin.next().strip()
        array = line.split()
        i = int(array[0])
        j = int(array[1])
        lens[i] = j
except:
    info = sys.exc_info()
    print info[0],info[1]
    cin.close()

cout = open(prefix+"/Cyto_conclusion","w")
for i,j in links:
    cout.write(str(i)+"("+str(lens[i])+")\t"+str(j)+"("+str(lens[j])+")\t"+str(flag[str(i)])+str(flag[str(j)])+"\n")
    
    le = len(emdis[(i,j)])
    cout.write(" emdis: -"+str(le)+" ")
    if le == 1:
        cout.write(str(emdis[(i,j)][0])+"\n")
    else:
        cout.write(str(emdis[(i,j)][0])+" / "+str(emdis[(i,j)][1])+"\n")
    
    le = len(refdis[(i,j)])
    cout.write("refdis: -"+str(le)+" ")
    if le == 1:
        cout.write(str(refdis[(i,j)][0]))
    else:
        cout.write(str(refdis[(i,j)][0])+" / "+str(refdis[(i,j)][1]))
    cout.write("\n")
    
    cout.write("Gap: ")
    if le == 1:
        cout.write(str(float(refdis[(i,j)][0]-lens[i])))
    else:
        cout.write(str(float(refdis[(i,j)][0]-lens[i]))+" / "+str(float(refdis[(i,j)][1]-lens[i])))
    cout.write("\n")

    test = 0
    f = 0
    cout.write("Position_i: ")
    for p in ctg_pos[str(i)]:
        cout.write(str(p)+" ")
        if abs(test-p) <= 1000:
            f = 1
        else:
            test = p    
    cout.write("F"+str(f)+"\n")
    f = 0
    test = 0
    cout.write("Position_j: ")
    for p in ctg_pos[str(j)]:
        cout.write(str(p)+" ")
        if abs(test-p) <= 1000:
            f = 1
        else:
            test = p    
    cout.write("F"+str(f)+"\n")

    cout.write("links dis("+str(len(links[(i,j)]))+"): ")
    for p in links[(i,j)]:
        cout.write(str(p)+" ")
    cout.write("\n")
    cout.write("\n")

cout.close()
