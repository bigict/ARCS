#! /usr/bin/python

Usage = '''
Author : Zheng QuanGang
Function : 

$python merging.py -I <insert size> 
		   -K <kmer> 
		   -C <cdbg> 
		   #-M <component> 
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

#K = K
line = ''
mer = ''
cov = []
lens = []
deltmp = []
#===================================================
# Reading cdbg file into contigs , seprate R and U , building kmer hash for trainning insert-size
#===================================================
cin = open(cdbg,"r")
try:
	while 1:
		line = cin.next().strip()
		array = line.split()
                cov.append(int(array[1]))
		line = cin.next().strip()
		contigs.append(line)
                lens.append(len(line))
		if len(line) > 2*I:
			for i in range(len(line)-K+1):
				tmp = line[i:i+K]
				if tmp in kmers:
					deltmp.append(tmp)
				else:
					kmers[tmp] = (index,i)
		index += 1
except:
	info = sys.exc_info()
	print info[0],info[1]
	cin.close()
	print "Number of kmers = " + str(len(kmers))
	print "[Done] Reading cdbg file into contigs , seprate R and U"
for i in deltmp:
        del kmers[tmp]
del deltmp
'''
	Reading reads file and trainning insert-size
'''
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
		if count > 1000:
			break
		temp1 = cin1.next().strip()
		temp1 = cin1.next().strip()
		temp2 = cin2.next().strip()
		temp2 = cin2.next().strip()
except:
    cin1.close()
    cin2.close()

em3.em3(isize)
print em3.mu1,em3.mu2
if em3.mu2 >= em3.mu1:
	Ireal = em3.mu2
	Idelta = em3.delta2
else:
	Ireal = em3.mu1
	Idelta = em3.delta1
II = int(2*int(Ireal) + 3*int(Idelta))
'''
II = 1000
Ireal = 179
Idelta = 10.5
print "II",II
print "Ireal",Ireal
print "Idelta",Idelta
print "[Done] Reading reads file and trainning insert-size"

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
		
del kmers
#del isize

#===================================================
#	Building kmer hash of U,R
#===================================================
delkmerin = {}
delkmerout = {}
kmerin = {}
kmerout = {}
index = 0
le = 0
tmp = 0
repeats = []
skip = 0
for i in range(len(contigs)):
        #if lens[i] < 80:
        #        continue
        if contigs[i][0:K] in refpos:
            if len(refpos[contigs[i][0:K]]) > 1:
                    R.append(i)
                    repeats.append(int(len(contigs[i])))
            else:
                    U.append(i)
        else:
            skip += 1
            continue
        '''
        if cov[i] == 1:
            if contigs[i][0:K] in refpos:
                if len(refpos[contigs[i][0:K]]) > 1:
                        R.append(i)
                else:
                        U.append(i)
            else:
                continue
                #tmp += 1
                #if tmp > 0:
                #    continue
                #U.append(i)
        elif cov[i] > 1:
            if contigs[i][0:K] in refpos:
                if len(refpos[contigs[i][0:K]]) > 1:
                        R.append(i)
                else:
                        U.append(i)
            else:
                continue
        else:
                continue
	'''
        temp = contigs[i]
	minl = min(II,len(temp))
        minl = len(temp)
        assert minl == lens[i]
        #print str(minl)
	for j in range(minl-K+1):
		if temp[j:j+K] in kmerin:
			delkmerin[temp[j:j+K]] += 1
		else:
			kmerin[temp[j:j+K]] = (i,j)
			delkmerin[temp[j:j+K]] = 1
	for j in range(len(temp)-minl,len(temp)-K+1):
		if temp[j:j+K] in kmerout:
			delkmerout[temp[j:j+K]] += 1
		else:
			kmerout[temp[j:j+K]] = (i,j)	
			delkmerout[temp[j:j+K]] = 1
	#for j in range(-minl,-K):
	#	if temp[j:j+K] in kmerout:
	#		delkmerout[temp[j:j+K]] += 1
	#	else:
	#		kmerout[temp[j:j+K]] = (i,len(temp)+j)	
	#		delkmerout[temp[j:j+K]] = 1
print "Number of R = " + str(len(R))
print "Number of U = " + str(len(U))
print "Number of all input k-mer hash = " + str(len(kmerin))
print "Number of all output k-mer hash = " + str(len(kmerout))
print "[Done] Building kmer hash of R"	
print "\t contigs not found in ref = "+str(skip)

testi = 0
try:
	for i in delkmerin:
		if int(delkmerin[i]) > 1:
			print "D " + i
			del kmerin[i]
			testi += 1
	for i in delkmerout:
		if int(delkmerout[i]) > 1:
			print "D " + i
			del kmerout[i]
			testi += 1
except:
	info = sys.exc_info()
	print info[0],info[1]
print "Deleting kmers (cov > 1) =" + str(testi)
print "Number of input k-mer hash = " + str(len(kmerin))
print "Number of output k-mer hash = " + str(len(kmerout))
print "[Done] Building kmer hash"	

repeats.sort()
count = 0
cout = open(prefix+"/Cyto_repeat_lens","w")
for i in repeats:
    cout.write(str(i)+"\n")
    if i >=100:
        count += 1
cout.close()
print "\t Number of Repeat = "+str(len(repeats))
print "\tRepeat length >= 100bp = "+str(count)
del repeats

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
count = 0
allc = 0
try:
	while 1:
		temp1 = cin1.next().strip()
		temp1 = cin1.next().strip()
		temp2 = cin2.next().strip()
		temp2 = cin2.next().strip()
		R2 = reversecomplement(temp2)
		minl = min(len(temp1),len(temp2))
                assert minl == 36
		for i in range(minl-K+1):
			kmer1 = temp1[i:i+K]
			kmer2 = R2[i:i+K]
			if kmer1 in kmerout and kmer2 in kmerin:
                                k1 = kmerout[kmer1]
                                k2 = kmerin[kmer2]
				if k1[0] != k2[0]:
                                        if (k1[0],k2[0]) in links:
						links[(k1[0],k2[0])].append(k1[1]+Ireal-k2[1])
						preads[(k1[0],k2[0])] += 1
						#preads[(kmerout[kmer1][0],kmerin[kmer2][0])].append(1) #(1,1) paired reads, paired k-mer
						#preads[(kmerout[kmer1][0],kmerin[kmer2][0])].append(1)
					else:
						links[(k1[0],k2[0])] = []
						links[(k1[0],k2[0])].append(k1[1]+Ireal-k2[1])
						#if flag == 1:
						preads[(k1[0],k2[0])] = 1
						#else:
						#	flag = 0
						#	preads[(kmerout[kmer1][0],kmerin[kmer2][0])][0] += 1
						#	preads[(kmerout[kmer1][0],kmerin[kmer2][0])][1] += 1
		R1 = reversecomplement(temp1)
		for i in range(minl-K+1):
			kmer1 = temp2[i:i+K]
			kmer2 = R1[i:i+K]
			if kmer1 in kmerout and kmer2 in kmerin:
                                k1 = kmerout[kmer1]
                                k2 = kmerin[kmer2]
				if k1[0] != k2[0]:
                                        if (k1[0],k2[0]) in links:
						links[(k1[0],k2[0])].append(k1[1]+Ireal-k2[1])
						preads[(k1[0],k2[0])] += 1
						#preads[(kmerout[kmer1][0],kmerin[kmer2][0])].append(1) #(1,1) paired reads, paired k-mer
						#preads[(kmerout[kmer1][0],kmerin[kmer2][0])].append(1)
					else:
						links[(k1[0],k2[0])] = []
						links[(k1[0],k2[0])].append(k1[1]+Ireal-k2[1])
						#if flag == 1:
						preads[(k1[0],k2[0])] = 1
						#else:
						#	flag = 0
						#	preads[(kmerout[kmer1][0],kmerin[kmer2][0])][0] += 1
						#	preads[(kmerout[kmer1][0],kmerin[kmer2][0])][1] += 1
		temp1 = cin1.next().strip()
		temp1 = cin1.next().strip()
		temp2 = cin2.next().strip()
		temp2 = cin2.next().strip()
		flag = 0
                count += 1
                print "\r"+str(count),
except:
	info = sys.exc_info()
	print info[0],info[1]
    	cin1.close()
    	cin2.close()
print "\n"
print "Number of links = " + str(len(links))
print "Number of reads = " + str(count)
print "[Done] Mapping paried reads to construct links"	


#===================================================
# filter : paired reads <= 1 and paired k-mer <= 40 
#===================================================
temp = []
for i in preads:
	#if preads[i][0] <= 0 or preads[i][1] <= 40:
        assert len(links[i]) == preads[i]
	if preads[i] <= 1000 or preads[i] >= 1400:
		#del preads[i]
		temp.append(i)
		del links[i]
for i in temp:
	del preads[i]
print "Number of deleted links = " + str(len(temp))
del temp
print "[Done] Delete links of paired reads < 0 or paired k-mer <= 40"

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
#	output paired-kmer dis
#===================================================
cout = open(prefix+"/Cyto_links_dis","w")

cout.write("i\tlen(i)\tj\tlen(j)\tRR/UR\t(num)\tdis...\n")
for i,j in links:
    if i in U and j in U:
        continue
    f1 = 'U'
    f2 = 'U'
    if i in R:
        f1 = 'R'
    if j in R:
        f2 = 'R'
    cout.write(str(i)+"\t"+str(len(contigs[i]))+"\t"+str(j)+"\t"+str(len(contigs[j]))+"\t"+f1+f2+"\t("+str(len(links[(i,j)]))+")\t")
    for k in links[(i,j)]:
        cout.write(str(k)+",")
    cout.write("\n")
cout.close()

#===================================================
#	estimate distance
#===================================================
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

emdis = {}
refdis = {}
diff = []
td = {}
thred = 50
linksnum = 0
correctnum = 0
uu = 0
uuc = 0
ur = 0
urc = 0
rr = 0
rrc = 0
count = 1
cout0 = open(prefix+"/Cyto_emdis","w")
cout1 = open(prefix+"/Cyto_refdis","w")
for i,j in links:
        print "count="+str(count)
        if i in U and j in U:
                #print "\r"+str(count)+"/"+str(len(links)),
                count += 1
                continue
        print "\tbefore em3"
	em3.em3(links[(i,j)])
        print "\tafter em3"
	mu1 = em3.mu1
	mu2 = em3.mu2
        #print "\r"+str(count)+"/"+str(len(links)),
        count += 1
	#print mu1,mu2
        #if mu1 > 0 and mu2 > 0 and abs(mu1 - mu2) <= 50:
        #    mu1 = (mu1+mu2)/2
        #    mu2 = 0
        if int(mu2) > 0 and int(mu1) == 0:
            mu1 = mu2
            mu2 = 0
        min1 = 0
	min2 = 0
	diff1 = 0
	diff2 = 0
	emdis[(i,j)] = (mu1,mu2)
        #start1 = (len(contigs[i])-K)/2
	start1 = 0
        kmer1 = contigs[i][start1:start1+K]
        #start2 = (len(contigs[j])-K)/2
	start2 = 0
        kmer2 = contigs[j][start2:start2+K]
	td[(i,j)] = []
	if kmer1 in refpos and kmer2 in refpos:
		for pi in refpos[kmer1]:
			for pj in refpos[kmer2]:
                                if pi >= 0 and pj >= 0 or pi < 0 and pj < 0:
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
				correctnum = 1
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
				correctnum = 1
		lable1 = "U"
                lable2 = "U"
                if i in U:
                        label1 = "U"
                        if j in U:
                                lable2 = "U"
                                uu += 1
                                uuc += correctnum
                        else:
                                lable2 = "R"
                                ur += 1
                                urc += correctnum 
                else:
                        label1 = "R"
                        if j in U:
                                lable2 = "U"
                                ur += 1
                                urc += correctnum
                        else:
                                lable2 = "R"
                                rr += 1
                                rrc += correctnum 
                correctnum = 0
		#diff[(i,j)] = (diff1,diff2)
		refdis[(i,j)] = (min1,min2)
		if (i,j) in refdis:
                        cout0.write(str(i)+"\t"+str(lable1)+"\t"+str(j)+"\t"+str(lable2)+"\t"+str(mu1)+","+str(mu2)+"\n")
			cout1.write(str(i)+"\t"+str(lable1)+"\t"+str(j)+"\t"+str(lable2)+"\t"+str(refdis[(i,j)][0])+","+str(refdis[(i,j)][1])+"\n")
cout0.close()
cout1.close()
#===================================================
# output Cyto information
#===================================================
'''
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
	if i in R:
		lable1 = "R"
	else:
		lable1 = "U"
	if j in R:
		lable2 = "R"
	else:
		lable2 = "U"

	if (i,j) in refdis:
		counta += 1
		if emdis[(i,j)][0] > 0 or emdis[(i,j)][1] > 0:
			if len(td[(i,j)]) > 1:
                                realD += 1
                                #for li in td[(i,j)]:
                                        #if int(li) < 2000: 
                                         #       print i,j
                                         #       break
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
print "[TEST] Correct dis of UU (<= 50)= " + str(uuc) + "/\t" +str(uu)
print "[TEST] Correct dis of UR (<= 50)= " + str(urc) + "/\t" +str(ur)
print "[TEST] Correct dis of RR (<= 50)= " + str(rrc) + "/\t" +str(rr)

print "[TEST] The number of correct links that truly have more than on dis = " + str(realD) + "/" + str(countd)
print "[TEST] The variance of diffs between emdis and refdis = " + str(mu)
print "[TEST] The standard deviation of diffs between emdis and refdis = " + str(math.sqrt(delta))
cout0.close()
cout1.close()
cout2.close()
'''	
print "[Done] estimate distances."

#===================================================
# print links for test
#===================================================
cout = open(prefix+"/Cyto_merging","w")
cout1 = open(prefix+"/Cyto_node","w")
for i,j in links:
        if i in U and j in U:
                continue
	if i in R:
		lable1 = "R"
	else:
		lable1 = "U"
	len1 = lens[i]
	if j in R:
		lable2 = "R"
	else:
		lable2 = "U"
	len2 = lens[j]
	#cout.write(str(i)+"|"+str(len1)+"\t"+str(lable1)+"\t"+str(j)+"|"+str(len2)+"\t"+str(lable2)+"\t"+str(preads[(i,j)][0])+"|"+str(preads[(i,j)][1])+"|"+str(dis[(i,j)])+"\n")
	cout.write(str(i)+"\t"+str(lable1)+"\t"+str(j)+"\t"+str(lable2)+"\t"+str(preads[(i,j)])+"|"+str(preads[(i,j)])+"|("+str(emdis[(i,j)][0])+","+str(emdis[(i,j)][1])+")\n")
cout.close()
index = 0
for i in range(len(contigs)):
	if i in R:
		label = "R"
	else:
		label = "U"
	cout1.write(str(i)+"\t"+str(lens[i])+"\t"+str(label)+"\n")
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






