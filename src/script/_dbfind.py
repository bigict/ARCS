#! /usr/bin/python

Usage = '''
Author : Zheng QuanGang
Function : 

$python dbfind.py   -C <condensed de Bruijn graph> 
                    -B <cdbg>
		    -D <workspace>
                    <start contig> <end contig>
'''

import sys,getopt
import os,math
import math

LEN = 50
opts, args = getopt.getopt(sys.argv[1:],"B:K:C:S:D:M:R:")
if len(opts) < 0:
	print Usage
	sys.exit()
if len(args) != 2:
	print Usage
	sys.exit()
else:
	s1 = int(args[0])
	s2 = int(args[1])
	
for op, value in opts:
	if op == '-K':
		K = int(value)
	elif op == '-C':
		debruijn = value
	elif op == '-S':
		scf = value
	elif op == '-D':
		prefix = value
	elif op == '-M':
		com = value
	elif op == '-R':
		ref = value
	elif op == '-B':
		cdbg_file = value
	else:
		print Usage
		sys.exit()

'''
	Function : Revese sequence
'''
def reversecomplement( read ):
	rcread = ''
	for i in range(len(read)):
		if read[len(read)- i - 1] == 'A':
			rcread = rcread + 'T'
			continue
		if read[len(read)- i - 1] == 'T':
			rcread = rcread + 'A'
			continue
		if read[len(read)- i - 1] == 'C':
			rcread = rcread + 'G'
			continue
		if read[len(read) - i - 1] == 'G':
			rcread = rcread + 'C'
			continue
		if read[len(read) - i - 1] == 'N':
			rcread = rcread + 'N'
        return rcread

#node_hash = [
de = {}
rcontigs = []
contigs = []
ctg_cov = []
ref_pos = {}
ctg_pos = []

cin = open(cdbg_file,"r")
try:
    while 1:
        line =cin.next().strip()
        line = cin.next().strip()
        rcontigs.append(line)
except:
    cin.close()
#===================================================
#	reading de Bruijn graph
#===================================================
l1 = 0 
l2 = 0
l3 = 0
R = []
U = []
index = 0
cindex = 0 
cin = open(debruijn,"r")
try:
    while 1:
        line = cin.next().strip()
        array = line.split()
        i = int(array[0])
        j = int(array[1])
        contigs.append(array[2])
        if array[2] in rcontigs:
            de[(i,j)] = index
            index += 1
            cindex += 1
        else:
            de[(i,j)] = "x"+str(cindex)
            cindex += 1

        line = cin.next().strip()
        array = line.split()
        ctg_cov.append(array)
        
except:
    info = sys.exc_info()
    print info[0],info[1]
    cin.close()
print "Step 1\nNumber of contigs = "+str(index)

index = 0
node_index = {}
for i,j in de:
    if i in node_index:
        node_index[i].append(j)
    else:
        node_index[i] = []
        node_index[i].append(j)
'''
de_nindex = {}
de_pindex = {}
for i,j in de:
    index = de[(i,j)]
    if index not in de_nindex:
        de_nindex[index] = []
    if j in node_index:
        for k in node_index[j]:
            de_nindex[index].append(de[(j,k)])
    index += 1
'''
de1node = {}
de2node = {}
for i,j in de:
    if type(de[(i,j)]) == int:
        de1node[int(de[(i,j)])] = i
        de2node[int(de[(i,j)])] = j
#===================================================
# building kmer hash for trainning insert-size
#===================================================
cin = open(prefix+"/Cyto_s",'w')
stack = []
start = de2node[s1]
end = de1node[s2]

print de1node[s1],de2node[s1],s1,len(rcontigs[s1])
print de1node[s2],de2node[s2],s2,len(rcontigs[s2])
cin.write("i\tj\tcontig_id|len\n")
cin.write(str(de1node[s1])+"\t"+str(de2node[s1])+"\t"+str(s1)+"|"+str(len(rcontigs[s1]))+"\n")
cin.write(str(de1node[s2])+"\t"+str(de2node[s2])+"\t"+str(s2)+"|"+str(len(rcontigs[s2]))+"\n")
stack.append(start)
depth = 0
count = 1
tmp = 0
while len(stack) != 0:
    k = stack.pop(0)
    count -= 1
    if k == end:
        print "found"
        if count <= 0:
            count = tmp
            tmp = 0
            depth += 1
        continue
    if depth >= 30:
        continue
    if k in node_index:
        for z in node_index[k]:
            if type(de[(k,z)]) == int:
                cin.write(str(k)+"\t"+str(z)+"\t"+str(de[(k,z)])+"|"+str(len(rcontigs[de[(k,z)]]))+"\n")
            else:
                index = int(de[(k,z)][1:])
                cin.write(str(k)+"\t"+str(z)+"\t"+str(de[(k,z)])+"|"+str(len(contigs[index]))+"\n")

            stack.append(z)
            tmp += 1
    if count <= 0:
        count = tmp
        tmp = 0
        depth += 1
    
cin.close()
    
    
    


