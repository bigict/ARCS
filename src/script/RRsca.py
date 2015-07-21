#! /usr/bin/python

Usage = '''
Author : Zheng QuanGang
Function : using repeats and scaffolds to construct graph

$python RRsca.py -I <insert size> 
		   -K <kmer> 
		   -C <cdbg> 
		   -S <scaffold> 
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

contigs = []
index = 0
kmers = {}
R = []
U = []

k1 = K
line = ''
mer = ''

#===================================================
# Reading cdbg file into contigs , seprate R and U , building kmer hash for trainning insert-size
#===================================================
cin = open(cdbg,"r")
try:
	while 1:
		line = cin.next().strip()
		array = line.split()
		if int(array[1]) > 1:
			R.append(index)
		else:
			U.append(index)
		line = cin.next().strip()
		contigs.append(line)
		if len(line) > 2*I:
			for i in range(len(line)-k1+1):
				tmp = line[i:i+K]
				if tmp in kmers:
					del kmers[tmp]
				else:
					kmers[tmp] = (index,i)
		index += 1
except:
	info = sys.exc_info()
	print info[0],info[1]
	cin.close()
	print "Number of R = " + str(len(R))
	print "Number of U = " + str(len(U))
	print "Number of kmers = " + str(len(kmers))
	print "[Done] Reading cdbg file into contigs , seprate R and U"

'''
	Reading reads file and trainning insert-size
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
				if abs(kmers[kmer2][1] - kmers[kmer1][1]) > I/5 and abs(kmers[kmer2][1] - kmers[kmer1][1]) < 2*I:
					isize.append(abs(kmers[kmer2][1] - kmers[kmer1][1]))
					count += 1
					break
		R1 = reversecomplement(temp1)
		for i in range(minl-K+1):
			kmer1 = temp2[i:i+K]
			kmer2 = R1[i:i+K]
			if kmer1 in kmers and kmer2 in kmers and kmers[kmer1][0] == kmers[kmer2][0]:
				if abs(kmers[kmer2][1] - kmers[kmer1][1]) > I/5 and abs(kmers[kmer2][1] - kmers[kmer1][1]) < 2*I:
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
    cin1.close()
    cin2.close()

em3.em3(isize)
if em3.mu2 != 0:
	Ireal = em3.mu2
	Idelta = em3.delta2
else:
	Ireal = em3.mu1
	Idelta = em3.delta1
II = int(Ireal + 3*Idelta)
II = 500
print Ireal
print Idelta
print "[Done] Reading reads file and trainning insert-size"

del kmers
del isize


delkmerin = {}
delkmerout = {}
component = []
lens = []  
kmerin = {}
kmerout = {}
index = 0
le = 0
#===================================================
#	Building kmer hash of R
#===================================================
'''
for i in R:
	lens.append(len(contigs[i]))
	component.append(('R',i,len(contigs[i])))
	temp = contigs[i]
	minl = min(II,len(temp))
	for j in range(minl-k1+1):
		if temp[j:j+k1] in kmerin:
			delkmerin[temp[j:j+k1]] += 1
		else:
			kmerin[temp[j:j+k1]] = (index,j,i)
			delkmerin[temp[j:j+k1]] = 1
	for j in range(-minl,-k1):
		if temp[j:j+k1] in kmerout:
			delkmerout[temp[j:j+k1]] += 1
		else:
			kmerout[temp[j:j+k1]] = (index,len(temp)+j,i)	
			delkmerout[temp[j:j+k1]] = 1
	index += 1	
print "Number of R = " + str(index)
print "Number of R input k-mer hash = " + str(len(kmerin))
print "Number of R output k-mer hash = " + str(len(kmerout))
print "[Done] Building kmer hash of R"	
'''
#===================================================
#	Reading scaffold and building kmer hash of U
#===================================================
cin= open(scf,"r")
scaff = []
count = 0
try:
	while 1:
		line = cin.next().strip()
		line = cin.next().strip()
                scaff.append(line)
		le = len(line)
		lens.append(le)
		component.append(("U",count,le))
                count += 1
		minl = min(II,len(line))
		for i in range(minl-k1+1):
			if line[i:i+k1] in kmerin:
				 delkmerin[line[i:i+k1]] += 1
			else:
				 kmerin[line[i:i+k1]] = (index,i,count)
				 delkmerin[line[i:i+k1]] = 1
		for i in range(-minl,-k1):
			if line[i:i+k1] in kmerout:
				delkmerout[line[i:i+k1]] += 1
			else:
				kmerout[line[i:i+k1]] = (index,le+i,count)			
				delkmerout[line[i:i+k1]] = 1
		index += 1
except:
	info = sys.exc_info()
	print info[0],info[1]
	cin.close()

print "NUmber of scaffolds = "+ str(count)
print "Number of input k-mer hash = " + str(len(kmerin))
print "Number of output k-mer hash = " + str(len(kmerout))
print "[Done] Reading scaffold and building kmer hash of all"
########  test ...
'''
cout = open(prefix+"/lens","w")
index = 0
for i in lens:
	cout.write(str(index)+"\t"+str(i)+"\n")
	index += 1
cout.close()
cout = open(prefix+"/kmerout","w")
for i in kmerout:
	cout.write(str(i)+"\t"+str(kmerout[i])+"\t"+str(delkmerout[i])+"\n")
cout.close()
cout = open(prefix+"/kmerin","w")
for i in kmerin:
	cout.write(str(i)+"\t"+str(kmerin[i])+"\t"+str(delkmerin[i])+"\n")
cout.close()
'''


testi = 0
try:
	for i in delkmerin:
		if int(delkmerin[i]) > 1:
			#print "D " + str(len(i))  + i
			del kmerin[i]
			testi += 1
	for i in delkmerout:
		if int(delkmerout[i]) > 1:
			#print "D " + str(len(i)) + i
			del kmerout[i]
			testi += 1
except:
	info = sys.exc_info()
	print info[0],info[1]
print "Deleting kmers (cov > 1) =" + str(testi)
print "Number of input k-mer hash = " + str(len(kmerin))
print "Number of output k-mer hash = " + str(len(kmerout))
print "[Done] Building kmer hash"	

del delkmerin
del delkmerout

#===================================================
#	 Mapping paried reads to construct links
#===================================================
testi = 0
testj = 0
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
				if kmerout[kmer1][0] != kmerin[kmer2][0]:
					if (kmerout[kmer1][0],kmerin[kmer2][0]) not in links:
						links[(kmerout[kmer1][0],kmerin[kmer2][0])] = []
						links[(kmerout[kmer1][0],kmerin[kmer2][0])].append(kmerout[kmer1][1]+Ireal-kmerin[kmer2][1])
						preads[(kmerout[kmer1][0],kmerin[kmer2][0])] = []
						preads[(kmerout[kmer1][0],kmerin[kmer2][0])].append(1) #(1,1) paired reads, paired k-mer
						preads[(kmerout[kmer1][0],kmerin[kmer2][0])].append(1)
					else:
						links[(kmerout[kmer1][0],kmerin[kmer2][0])].append(kmerout[kmer1][1]+Ireal-kmerin[kmer2][1])
						if flag == 1:
							preads[(kmerout[kmer1][0],kmerin[kmer2][0])][1] += 1
						else:
							flag = 1
							preads[(kmerout[kmer1][0],kmerin[kmer2][0])][0] += 1
							preads[(kmerout[kmer1][0],kmerin[kmer2][0])][1] += 1
				
		R1 = reversecomplement(temp1)
		for i in range(minl-K+1):
			kmer1 = temp2[i:i+K]
			kmer2 = R1[i:i+K]
			if kmer1 in kmerout and kmer2 in kmerin:
				if kmerout[kmer1][0] != kmerin[kmer2][0]:
					if (kmerout[kmer1][0],kmerin[kmer2][0]) not in links:
						links[(kmerout[kmer1][0],kmerin[kmer2][0])] = []
						links[(kmerout[kmer1][0],kmerin[kmer2][0])].append(kmerout[kmer1][1]+Ireal-kmerin[kmer2][1])
						preads[(kmerout[kmer1][0],kmerin[kmer2][0])] = []
						preads[(kmerout[kmer1][0],kmerin[kmer2][0])].append(1) #(1,1) paired reads, paired k-mer
						preads[(kmerout[kmer1][0],kmerin[kmer2][0])].append(1)
					else:
						links[(kmerout[kmer1][0],kmerin[kmer2][0])].append(kmerout[kmer1][1]+Ireal-kmerin[kmer2][1])
						if flag == 1:
							preads[(kmerout[kmer1][0],kmerin[kmer2][0])][1] += 1
						else:
							flag = 1
							preads[(kmerout[kmer1][0],kmerin[kmer2][0])][0] += 1
							preads[(kmerout[kmer1][0],kmerin[kmer2][0])][1] += 1
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
print "Number of links = " + str(len(links))
print "[Done] Mapping paried reads to construct links"	


#===================================================
# filter : paired reads <= 1 and paired k-mer <= 40 
#===================================================
temp = []
for i in preads:
	if preads[i][0] <= 1 or preads[i][1] <= 0:
		#del preads[i]
		temp.append(i)
		del links[i]
for i in temp:
	del preads[i]
print "Number of deleted links = " + str(len(temp))
del temp
print "[Done] Delete links of paired reads <= 1 or paired k-mer <= 0"

#===================================================
# building links index
#===================================================
lindex = {}
for i,j in links:
	if i in lindex:
		lindex[i].append(j)
	else:
		lindex[i] = []
		lindex[i].append(j)	
print "[Done] building links index"	

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
refpos = {}
for i in range(len(refs)-K+1):
	kmer = refs[i:i+K]
	if kmer in refpos:
		refpos[kmer].append(i)
	else:
		refpos[kmer] = []
		refpos[kmer].append(i)
for i in range(len(ref_r)-K+1):
	kmer = ref_r[i:i+K]
	if kmer in refpos:
		refpos[kmer].append(-i)
	else:
		refpos[kmer] = []
		refpos[kmer].append(-i)
		
#===================================================
#	output paired-kmer dis
#===================================================
cout = open(prefix+"/Cyto_links_dis","w")
for i,j in links:
    cout.write(str(i)+"\t"+str(j)+"\t")
    for k in links[(i,j)]:
        cout.write(str(k)+",")
    cout.write("\n")
cout.close()

#===================================================
#	estimate distance
#===================================================
cout = open(prefix+"/Cyto_contig_pos","w")
for i in range(len(component)):
        if component[i][0] == "R":
        	kmer = contigs[int(component[i][1])][0:K]
        else:
                kmer = scaff[int(component[i][1])][0:K]
        if kmer in refpos:
		cout.write(str(i)+"\t")
		for j in refpos[kmer]:
			cout.write(str(j)+",")
		cout.write("\n")
cout.close()

emdis = {}
refdis = {}
diff = []
td = {}
thred = 100
linksnum = 0
correctnum = 0
for i,j in links:
	em3.em3(links[(i,j)]) 
	mu1 = em3.mu1
	mu2 = em3.mu2
	#print mu1,mu2
        if mu1 > 0 and mu2 > 0 and abs(mu1 - mu2) <= 70:
            mu1 = (mu1+mu2)/2
            mu2 = 0
	min1 = 0
	min2 = 0
	diff1 = 0
	diff2 = 0
	emdis[(i,j)] = (mu1,mu2)
        if component[i][0] == "R":
                start1 = (len(contigs[int(component[i][1])])-K)/2
	        kmer1 = contigs[int(component[i][1])][start1:start1+K]
        else:
                start1 = (len(scaff[int(component[i][1])])-K)/2
	        kmer1 = scaff[int(component[i][1])][start1:start1+K]
        if component[j][0] == "R":
                start2 = (len(contigs[int(component[j][1])])-K)/2
	        kmer2 = contigs[int(component[j][1])][start2:start2+K]
        else:
                start2 = (len(scaff[int(component[j][1])])-K)/2
	        kmer2 = scaff[int(component[j][1])][start2:start2+K]
	td[(i,j)] = []
	if kmer1 in refpos and kmer2 in refpos:
		for pi in refpos[kmer1]:
			for pj in refpos[kmer2]:
                                if (pi >= 0 and pj >= 0) or (pi < 0 and pj < 0):
                			td[(i,j)].append(abs(abs(pi)-start1-abs(pj)+start2))
                if len(td[i,j]) == 0:
                        continue
                if mu1 > 0:
			linksnum += 1
			diff1 = sys.maxint
			min1 = sys.maxint
			for m1 in td[(i,j)]:
				if abs(m1 - mu1) < diff1:
					diff1 = abs(m1-mu1)
					min1 = m1			
			diff.append(diff1)
			if diff1 <= thred:
				correctnum += 1
		if mu2 > 0:
			linksnum += 1
			diff2 = sys.maxint
			min2 = sys.maxint
			for m2 in td[(i,j)]:
				if abs(m2 - mu2) < diff2:
					diff2 = abs(m2-mu2)
					min2 = m2
			diff.append(diff2)
			if diff2 <= thred:
				correctnum += 1
		
		refdis[(i,j)] = (min1,min2)
cout = open(prefix+"/Cyto_alldisinref","w")
for i,j in td:
    cout.write(str(i)+"\t"+str(j)+"\t")
    for d in td[(i,j)]:
        cout.write(str(d)+",")
    cout.write("\n")
cout.close()
#===================================================
# output Cyto information
#===================================================
cout0 = open(prefix+"/Cyto_emdis","w")
cout1 = open(prefix+"/Cyto_refdis","w")
cout2 = open(prefix+"/Cyto_diff","w")
countd = 0
counta = 0
total = 0
realD = 0 # readl doubel paths
delta = 0.0
mu = 0.0
for i,j in emdis:
	total += 1
        lable1 = component[i][0]
        lable2 = component[j][0]

	if (i,j) in refdis:
		counta += 1
		if emdis[(i,j)][0] > 0 and emdis[(i,j)][1] > 0:
			if len(td[(i,j)]) > 1:
				realD += 1
			countd += 1
			cout0.write(str(i)+"\t"+str(lable1)+"\t"+str(j)+"\t"+str(lable2)+"\t"+str(emdis[(i,j)][0])+","+str(emdis[(i,j)][1])+"\n")
			cout1.write(str(i)+"\t"+str(lable1)+"\t"+str(j)+"\t"+str(lable2)+"\t"+str(refdis[(i,j)][0])+","+str(refdis[(i,j)][1])+"\n")
		elif emdis[(i,j)][0] > 0:	
			cout0.write(str(i)+"\t"+str(lable1)+"\t"+str(j)+"\t"+str(lable2)+"\t"+str(emdis[(i,j)][0])+"\n")
			cout1.write(str(i)+"\t"+str(lable1)+"\t"+str(j)+"\t"+str(lable2)+"\t"+str(refdis[(i,j)][0])+"\n")
		elif emdis[(i,j)][1] > 0:	
			cout0.write(str(i)+"\t"+str(lable1)+"\t"+str(j)+"\t"+str(lable2)+"\t"+str(emdis[(i,j)][1])+"\n")
			cout1.write(str(i)+"\t"+str(lable1)+"\t"+str(j)+"\t"+str(lable2)+"\t"+str(refdis[(i,j)][1])+"\n")
mu = float(sum(diff)/len(diff))
for i in diff:
	cout2.write(str(i)+"\n")
	delta += int(i)*int(i)	


print "[TEST] The number of all dis = " + str(total)
print "[TEST] The number of links found in ref = " + str(counta)
print "[TEST] The number of links that have two dis = " + str(countd)
print "[TEST] The percentage of correct dis (<= 100)= " + str(correctnum) +"\t"+ str(linksnum) + "\t" + str(float(correctnum)/linksnum)
print "[TEST] The number of correct links that truly have more than on dis = " + str(realD) + "/" + str(countd)
print "[TEST] The variance of diffs between emdis and refdis = " + str(mu)
print "[TEST] The standard deviation of diffs between emdis and refdis = " + str(math.sqrt(delta))
cout0.close()
cout1.close()
cout2.close()
	
print "[Done] estimate distances."

#===================================================
#	output component position in ref
#===================================================
cout = open(prefix+"/Cyto_refpos","w")
for i in range(len(component)):
        if component[i][0] == "R":
        	kmer = contigs[int(component[i][1])][0:K]
        else:
                kmer = scaff[int(component[i][1])][0:K]
	if kmer in refpos:
		cout.write(str(i)+"\t")
		for j in refpos[kmer]:
			cout.write(str(j)+",")
		cout.write("\n")
cout.close() 

#===================================================
# print links for test
#===================================================
cout = open(prefix+"/Cyto_merging","w")
cout1 = open(prefix+"/Cyto_node","w")
for i,j in links:
        lable1 = component[i][0]
	len1 = lens[i]
        lable2 = component[j][0]
	len2 = lens[j]
	#cout.write(str(i)+"|"+str(len1)+"\t"+str(lable1)+"\t"+str(j)+"|"+str(len2)+"\t"+str(lable2)+"\t"+str(preads[(i,j)][0])+"|"+str(preads[(i,j)][1])+"|"+str(dis[(i,j)])+"\n")
	cout.write(str(i)+"\t"+str(lable1)+"\t"+str(j)+"\t"+str(lable2)+"\t"+str(preads[(i,j)][0])+"|"+str(preads[(i,j)][1])+"|("+str(emdis[(i,j)][0])+","+str(emdis[(i,j)][1])+")\n")
cout.close()
index = 0
for i in range(len(component)):
        lable = component[i][0]
	cout1.write(str(i)+"\t"+str(lens[i])+"\t"+str(lable)+"\n")
cout1.close()

'''
#
#	Output position_lp.math
#
cout = open(prefix+"/position_lp_10.math","w")
count = 0
sep = 0
for i in range(len(component)):
	cout.write("var x_"+str(i)+";\n")
for i,j in links:
	cout.write("var e_"+str(i)+"_"+str(j)+";\n")
	cout.write("var E_"+str(i)+"_"+str(j)+";\n")
	count += 1

cout.write("\n")
cout.write("minimize z: ")
for i,j in links:
	cout.write("E_"+str(i)+"_"+str(j)+" + ")
	sep += 1
	if sep%10 == 0:
		cout.write("\n")
cout.write("0;\n\n")
sep = 1
for i,j in links:
	cout.write("s.t. con"+str(sep)+" : x_"+str(j)+" - x_"+str(i)+" + e_"+str(i)+"_"+str(j)+" = "+str(dis[(i,j)])+";\n")
	sep += 1
for i,j in links:
	cout.write("s.t. con"+str(sep)+" : E_"+str(i)+"_"+str(j)+" + e_"+str(i)+"_"+str(j)+" >= "+str(0)+";\n")
	sep += 1
	cout.write("s.t. con"+str(sep)+" : E_"+str(i)+"_"+str(j)+" - e_"+str(i)+"_"+str(j)+" >= "+str(0)+";\n")
	sep += 1

cout.write("\n\nend;\n")
'''

#
# deal with R - R
#






