#! /usr/bin/python

Usage = '''
Author : Zheng QuanGang
Function : 

$python deBruijn.py -I <insert size> 
		    -K <k-mer> 
		    -C <condensed de Bruijn graph> 
		    -D <workspace>
		    -R <ref>
		    <reads 1> <reads 2>
'''

import sys,getopt
import os,math
import em3
import math

LEN = 50
opts, args = getopt.getopt(sys.argv[1:],"I:K:C:S:D:M:R:")
if len(opts) < 0:
	print Usage
	sys.exit()
if len(args) != 2:
	print Usage
	sys.exit()
else:
	rds1 = args[0]
	rds2 = args[1]
	
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
	elif op == '-I':
		I = int(value)
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
contigs = []
ctg_cov = []
ref_pos = {}
ctg_pos = []

#===================================================
#	ref to kmer hash
#===================================================
ref_in = open(ref, "r")
refs = ''
ref_in.next()
try:
	while 1:
		refs += ref_in.next().strip()
except:
	ref_in.close()
ref_r = reversecomplement(refs)
for i in range(len(refs)-K+1):
	kmer = refs[i:i+K]
	if kmer in ref_pos:
		ref_pos[kmer].append(i)
	else:
		ref_pos[kmer] = []
		ref_pos[kmer].append(i)
for i in range(len(ref_r)-K+1):
	kmer = ref_r[i:i+K]
	if kmer in ref_pos:
		ref_pos[kmer].append(-i)
	else:
		ref_pos[kmer] = []
		ref_pos[kmer].append(-i)
#===================================================
#	reading de Bruijn graph
#===================================================
l1 = 0 
l2 = 0
l3 = 0
R = []
U = []
index = 0
cin = open(debruijn,"r")
try:
    while 1:
        line = cin.next().strip()
        array = line.split()
        i = int(array[0])
        j = int(array[1])
        de[(i,j)] = index
        contigs.append(array[2])
        kmer = array[2][0:K]
        tmp = []
        if kmer in ref_pos:
            tmp = ref_pos[kmer]
        ctg_pos.append(tmp)
        if len(tmp) == 1:
            U.append(index)
        elif len(tmp) > 1:
            R.append(index)

        if len(array[2]) <= 50:
            l1 += 1
        elif len(array[2]) <= 100:
            l2 += 1
        else:
            l3 += 1

        line = cin.next().strip()
        array = line.split()
        ctg_cov.append(array)
        
        index += 1
except:
    info = sys.exc_info()
    print info[0],info[1]
    cin.close()
print "Step 1\nNumber of contigs = "+str(index)
print "\t         <= 50bp = "+str(l1)
print "\t>50bp && <=100bp = "+str(l2)
print "\t         > 100bp = "+str(l3)
print "Number of unique contigs = "+str(len(U))
print "Number of repeat contigs = "+str(len(R))

index = 0
node_index = {}
for i,j in de:
    if i in node_index:
        node_index[i].append(j)
    else:
        node_index[i] = []
        node_index[i].append(j)
de_nindex = {}
de_pindex = {}
for i,j in de:
    if index not in de_nindex:
        de_nindex[index] = []
    if j in node_index:
        for k in node_index[j]:
            de_nindex[index].append(de[(j,k)])
    index += 1
#===================================================
# building kmer hash for trainning insert-size
#===================================================
print "Step 2\nTrainning insert size..."
kmers = {}
deltmp = {}
index = 0
for ctg in contigs:
    if len(ctg) > 2*I:
	for i in range(len(ctg)-K+1):
	    tmp = ctg[i:i+K]
	    if tmp in kmers:
	    	deltmp[tmp] = 1
	    else:
	    	kmers[tmp] = (index,i)
    index += 1

#print "Number of kmers = " + str(len(kmers))
for i in deltmp:
        del kmers[tmp]
del deltmp
'''	Reading reads file and trainning insert-size
'''
cin1 = open(rds1,"r")
cin2 = open(rds2,"r")
isize = []
count = 0
try:
	while 1:
		temp1 = cin1.next().strip()
		temp1 = cin1.next().strip()
		temp2 = cin2.next().strip()
		temp2 = cin2.next().strip()
		R2 = reversecomplement(temp2)
		minl = min(len(temp1),len(temp2))
		for i in range(minl-K+1):
			kmer1 = temp1[i:i+K]
			kmer2 = R2[i:i+K]
			if kmer1 in kmers and kmer2 in kmers and kmers[kmer1][0] == kmers[kmer2][0]:
				if abs(kmers[kmer2][1] - kmers[kmer1][1]) > I/2 and abs(kmers[kmer2][1] - kmers[kmer1][1]) < 2*I:
					isize.append(abs(kmers[kmer2][1] - kmers[kmer1][1]))
					count += 1
					break
		R1 = reversecomplement(temp1)
		for i in range(minl-K+1):
			kmer1 = temp2[i:i+K]
			kmer2 = R1[i:i+K]
			if kmer1 in kmers and kmer2 in kmers and kmers[kmer1][0] == kmers[kmer2][0]:
				if abs(kmers[kmer2][1] - kmers[kmer1][1]) > I/2 and abs(kmers[kmer2][1] - kmers[kmer1][1]) < 2*I:
					isize.append(abs(kmers[kmer2][1] - kmers[kmer1][1]))
					count += 1
					break 
		if count > 10000:
			break
		temp1 = cin1.next().strip()
		temp1 = cin1.next().strip()
		temp2 = cin2.next().strip()
		temp2 = cin2.next().strip()
except:
    info = sys.exc_info()
    print info[0],info[1]
    cin1.close()
    cin2.close()

em3.em3(isize)
print "\tEm = ",em3.mu1,em3.mu2
if em3.mu2 >= em3.mu1:
	Ireal = em3.mu2
	Idelta = em3.delta2
else:
	Ireal = em3.mu1
	Idelta = em3.delta1
II = int(2*Ireal + 3*Idelta)
print "\tinsert size = ",Ireal
print "\t     delata = ",Idelta
#print "[Done] Reading reads file and trainning insert-size"

#===================================================
#	Building kmer hash of all contigs
#===================================================
delkmerin = {}
delkmerout = {}
kmerin = {}
kmerout = {}
index = 0
le = 0
tmp = 0
for i in range(len(contigs)):
	temp = contigs[i]
	minl = min(II,len(temp))
	for j in range(minl-K+1):
		if temp[j:j+K] in kmerin:
			delkmerin[temp[j:j+K]] = 1
		else:
			kmerin[temp[j:j+K]] = (i,j)
	for j in range(-minl,-K):
		if temp[j:j+K] in kmerout:
			delkmerout[temp[j:j+K]] = 1
		else:
			kmerout[temp[j:j+K]] = (i,len(temp)+j)	
for i in delkmerin:
	print "D " + i
	del kmerin[i]
for i in delkmerout:
	print "D " + i
	del kmerout[i]
print "Step 3\nbuilding k-mer hash of all contigs"
print "\tNumber of input k-mer = " + str(len(kmerin))
print "\tNumber of output k-mer = " + str(len(kmerout))
del delkmerin
del delkmerout

#===================================================
#	 Mapping paried reads to construct links
#===================================================
links = {}
preads = {}

cin1 = open(rds1,"r")
cin2 = open(rds2,"r")
isize = []
flag = 0
try:
	while 1:
		temp1 = cin1.next().strip()
		temp1 = cin1.next().strip()
		temp2 = cin2.next().strip()
		temp2 = cin2.next().strip()
		R2 = reversecomplement(temp2)
		minl = min(len(temp1),len(temp2))
		for i in range(minl-K+1):
			kmer1 = temp1[i:i+K]
			kmer2 = R2[i:i+K]
			if kmer1 in kmerout and kmer2 in kmerin:
                                k1 = kmerout[kmer1]
                                k2 = kmerin[kmer2]
				if k1[0] != k2[0]:
					if (k1[0],k2[0]) not in links:
						links[(k1[0],k2[0])] = []
						links[(k1[0],k2[0])].append(k1[1]+Ireal-k2[1])
						preads[(k1[0],k2[0])] = []
						preads[(k1[0],k2[0])].append(1) #(1,1) paired reads, paired k-mer
						preads[(k1[0],k2[0])].append(1)
					else:
						links[(k1[0],k2[0])].append(k1[1]+Ireal-k2[1])
						if flag == 1:
							preads[(k1[0],k2[0])][1] += 1
						else:
							flag = 1
							preads[(k1[0],k2[0])][0] += 1
							preads[(k1[0],k2[0])][1] += 1
				
		R1 = reversecomplement(temp1)
		for i in range(minl-K+1):
			kmer1 = temp2[i:i+K]
			kmer2 = R1[i:i+K]
			if kmer1 in kmerout and kmer2 in kmerin:
                                k1 = kmerout[kmer1]
                                k2 = kmerin[kmer2]
				if k1[0] != k2[0]:
					if (k1[0],k2[0]) not in links:
						links[(k1[0],k2[0])] = []
						links[(k1[0],k2[0])].append(k1[1]+Ireal-k2[1])
						preads[(k1[0],k2[0])] = []
						preads[(k1[0],k2[0])].append(1) #(1,1) paired reads, paired k-mer
						preads[(k1[0],k2[0])].append(1)
					else:
						links[(k1[0],k2[0])].append(k1[1]+Ireal-k2[1])
						if flag == 1:
							preads[(k1[0],k2[0])][1] += 1
						else:
							flag = 1
							preads[(k1[0],k2[0])][0] += 1
							preads[(k1[0],k2[0])][1] += 1
				
		temp1 = cin1.next().strip()
		temp1 = cin1.next().strip()
		temp2 = cin2.next().strip()
		temp2 = cin2.next().strip()
		flag = 0
except:
	info = sys.exc_info()
	print info[0],info[1]
    	cin1.close()
    	cin2.close()
print "Step 4\n Mapping paired k-mer"

''' filter : paired reads &&  paired k-mer  
'''
temp = []
for i in preads:
	if preads[i][0] <= 2 and preads[i][1] <= 20:
		temp.append(i)
		del links[i]
for i in temp:
	del preads[i]
print "\tDelete links : paired reads < 2 and paired k-mer <= 20"
print "\tNumber of deleted links = " + str(len(temp))
del temp
print "\tNumber of links = " + str(len(links))

#===================================================
# building links index
#===================================================
nindex = {}
pindex = {}
for i,j in links:
	if i in nindex:
		nindex[i].append(j)
	else:
		nindex[i] = []
		nindex[i].append(j)	
	if j in pindex:
		pindex[j].append(i)
	else:
		pindex[j] = []
		pindex[j].append(i)	
print "\tbuilding links index"	
#===================================================
#  in-degree and out-degree
#===================================================
i = 0
in_degree = [0 for i in range(len(contigs))]
out_degree = [0 for i in range(len(contigs))]
for i,j in links:
	out_degree[i] += 1
	in_degree[j] += 1
print "\tbuilding in-degree and out-degree for each contig"	

#===================================================
#   estimate distance
#===================================================
'''
cout = open(prefix+"/Cyto_contig_pos","w")
for i in R:
	kmer = contigs[i][0:K]
	if kmer in refpos:
		cout.write(str(i)+"\tR\t")
		for j in refpos[kmer]:
			cout.write(str(j)+",")
		cout.write("\n")
for i in U:
	kmer = contigs[i][0:K]
	if kmer in refpos:
		cout.write(str(i)+"\tU\t")
		for j in refpos[kmer]:
			cout.write(str(j)+",")
		cout.write("\n")
cout.close()
'''
emdis = {}
refdis = {}
td = {}
for i,j in links:
	em3.em3(links[(i,j)])
	mu1 = int(em3.mu1)
	mu2 = int(em3.mu2)
        if mu2 > 0 and mu1 > 0:
            mu1 = (mu2 + mu1) / 2
        else:
            mu1 = mu1 + mu2
            mu2 = 0
        min1 = 0
	diff1 = 0
	emdis[(i,j)] = mu1
        #start1 = (len(contigs[i])-K)/2
	start1 = 0
        kmer1 = contigs[i][start1:start1+K]
        #start2 = (len(contigs[j])-K)/2
	start2 = 0
        kmer2 = contigs[j][start2:start2+K]
	td[(i,j)] = []
	if kmer1 in ref_pos and kmer2 in ref_pos:
		for pi in ref_pos[kmer1]:
			for pj in ref_pos[kmer2]:
                                if pi >= 0 and pj >= 0 or pi < 0 and pj < 0:
                			td[(i,j)].append(abs(abs(pi)-start1-abs(pj)+start2))
                if len(td[i,j]) == 0:
                        continue
                if mu1 > 0:
			diff1 = sys.maxint
			min1 = sys.maxint
			for m1 in td[(i,j)]:
				if abs(m1 - mu1) < diff1:
					diff1 = abs(m1-mu1)
					min1 = m1			
		
		refdis[(i,j)] = min1
print "Step 5\nestimate distance..."

#===================================================
#  split node in deBruijn based on links
#===================================================
dcount = 0
for i in de_nindex:
    for j in de_nindex[i]:
        if (i,j) not in links:
            del de_nindex[i][de_nindex[i].index(j)]
            dcount += 1

print "\tbreak connections in de Bruijn = "+str(dcount)
#===================================================
#  find paths
#===================================================
print "Step 6\n Using depth first search to find paths between i,j "
i = 0
paths = {}
directlink = 0
plens = {}
stack = []
dfsflag = [0 for i in contigs] #0:unvisited and unhandled 1:visited and unhandle 2:visited and handled 

def DFS(a,b,c,ltmp,ptmp,path,plen):
    if c == b:                     
        print "found"
        path.append(ptmp)          
        plen.append(ltmp)          
        return                    
    if ltmp+len(contigs[c]) >= emdis[(a,b)] + I:
        return
    ptmp.append(c)
    ltmp += len(contigs[c])
    for k in de_nindex[c]:             
        DFS(a,b,k,ltmp,ptmp,path,plen)
    tmp = ptmp.pop()
    ltmp -= len(contigs[c])


for i,j in links:
    if j in de_nindex[i]:
        directlink += 1
        break
    paths[(i,j)] = []
    plens[(i,j)] = []
    ptmp = []
    ltmp = 0#len(contigs[i])
    #ptmp.append(i)
    k = i
    DFS(i,j,k,ltmp,ptmp,paths[(i,j)],plens[(i,j)])
print "\tdirectlinks = ",directlink

#===================================================
# output Cyto information
#===================================================
print "Step 7\nOutputing Cyto information"
cout0 = open(prefix+"/Cyto_deBruijn","w")
cout1 = open(prefix+"/Cyto_scaf","w")
cout2 = open(prefix+"/Cyto_paths","w")
cout3 = open(prefix+"/Cyto_downsizing_deBruijn","w")
# i j index len
for i,j in de:
    cout0.write(str(i)+"\t"+str(j)+"\t"+str(de[(i,j)])+"\t"+str(len(contigs[de[(i,j)]]))+"\n")

cout3.write("I\tJ\tlen(I)|len(J)\n")
for i in de_nindex:
    for j in de_nindex[i]:
        cout3.write(str(i)+"\t"+str(j)+"\t"+str(len(contigs[i]))+"|"+str(len(contigs[j]))+"\n")

#i len j len RR/RU/UU em refdis dis
for i,j in links:
    if i in R:
        flag1 = 'R'
    else:
        flag1 = 'U'
    if j in R:
        flag2 = 'R'
    else:
        flag2 = 'U'
    cout1.write(str(i)+"\t"+str(len(contigs[i]))+"\t"+str(j)+"\t"+str(len(contigs[j]))+"\t"+str(flag1)+str(flag2)+"\t")
    for di in links[(i,j)]:
        cout1.write(str(di)+",")
    cout1.write("\n")

for i,j in paths:
    path = paths[(i,j)]
    for k in range(len(path)):
        cout2.write("%6d-%6d (%2d)\t"%(i,j,k))
        for e in path[k]:
            cout2.write(str(e)+",")
        cout2.write("\n")
cout0.close()
cout1.close()
cout2.close()
cout3.close()


