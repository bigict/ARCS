#! /usr/bin/python

Usage = '''
Author : Zheng QuanGang
Function : 

$python merging.py -I <insert size> 
		   -K <kmer> 
		   -C <cdbg> 
		   -S <scaffold> 
		   -M <component> 
		   -D <workspace>
		   <reads 1> <reads 2>
'''

import sys,getopt
import os,math
import em3

LEN = 50
opts, args = getopt.getopt(sys.argv[1:],"I:K:C:S:D:M:")
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

'''
 Reading cdbg file into contigs , seprate R and U , building kmer hash for trainning
 insert-size
'''
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
if em3.mu2 != 0:
	Ireal = em3.mu2
	Idelta = em3.delta2
else:
	Ireal = em3.mu1
	Idelta = em3.delta1
II = int(Ireal + 3*Idelta)
print Ireal
print Idelta
print "[Done] Reading reads file and trainning insert-size"

del kmers
del isize

#
#	Reading component and building kmer hash of U
#
delkmerin = {}
delkmerout = {}
component = []
lens = []  # component lens
gaps = []
kmerin = {}
kmerout = {}
cin= open(com,"r")
index = 0
le = 0
try:
	while 1:
		line = cin.next().strip()
		line = cin.next().strip()
		array = line.split()
		line = cin.next().strip()
		array1 = line.split()
		t = []
		g = []
		for i in array:
			t.append(i)
			le += len(contigs[int(i)]) - K + 1
		for i in array1:
			g.append(i)
			le += int(i)
		le += K
		lens.append(le)
		component.append((t,g))
		temp = contigs[int(array[0])]
		minl = min(II,len(temp))
		for i in range(minl-k1+1):
			if temp[i:i+k1] in kmerin:
				 delkmerin[temp[i:i+k1]] += 1
			else:
				 kmerin[temp[i:i+k1]] = (index,i,array[0])
				 delkmerin[temp[i:i+k1]] = 1
		temp = contigs[int(array[len(array)-1])]
		minl = min(II,len(temp))
		for i in range(-minl,-k1):
			if temp[i:i+k1] in kmerout:
				delkmerout[temp[i:i+k1]] += 1
			else:
				kmerout[temp[i:i+k1]] = (index,le+i,array[len(array)-1])			
				delkmerout[temp[i:i+k1]] = 1
		index += 1
		le = 0
except:
	info = sys.exc_info()
	print info[0],info[1]
	cin.close()

print "NUmber of scaffolds contained in component = "+ str(index)
print "Number of U input k-mer hash = " + str(len(kmerin))
print "Number of U output k-mer hash = " + str(len(kmerout))
print "[Done] Reading component and building kmer hash of U"

#
#	Building kmer hash of R
#
for i in R:
	t = []
	g = []
	t.append(i)
	lens.append(len(contigs[i]))
	component.append((t,g))
	temp = contigs[i]
	minl = min(II,len(temp))
	if len(temp) >= Ireal:   # threshold 2K or Ireal, II
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
print "Number of all input k-mer hash = " + str(len(kmerin))
print "Number of all output k-mer hash = " + str(len(kmerout))
print "[Done] Building kmer hash of R"	

########  test ...
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
print "[Done] Building kmer hash of R"	

del delkmerin
del delkmerout

#
#	 Mapping paried reads to construct links
#
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


#
# filter : paired reads <= 1 and paired k-mer <= 40 
#
temp = []
for i in preads:
	if preads[i][0] < 10 or preads[i][1] < 80:
		#del preads[i]
		temp.append(i)
		del links[i]
for i in temp:
	del preads[i]
print "Number of deleted links = " + str(len(temp))
del temp
print "[Done] Delete links of paired reads <= 1 and paired k-mer <= 40"

#
# building links index
#
lindex = {}
for i,j in links:
	if i in lindex:
		lindex[i].append(j)
	else:
		lindex[i] = []
		lindex[i].append(j)	
print "[Done] building links index"	

#
#	estimate distance
#
dis = {}
for i,j in links:
	d = int(sum(links[(i,j)])/len(links[(i,j)]))
	dis[(i,j)] = d
print "[Done] estimate distances."


# print links for test
cout = open(prefix+"/Cyto_merging","w")
cout1 = open(prefix+"/Cyto_node","w")
for i,j in links:
	if component[i][0][0] in R:
		lable1 = "R"
	else:
		lable1 = "U"
	len1 = lens[i]
	if component[j][0][0] in R:
		lable2 = "R"
	else:
		lable2 = "U"
	len2 = lens[j]
	cout.write(str(i)+"|"+str(len1)+"\t"+str(lable1)+"\t"+str(j)+"|"+str(len2)+"\t"+str(lable2)+"\t"+str(preads[(i,j)][0])+"|"+str(preads[(i,j)][1])+"|"+str(dis[(i,j)])+"\n")
cout.close()
index = 0
for i in range(len(component)):
	if component[i][0][0] in R:
		label = "R"
	else:
		label = "U"
	cout1.write(str(i)+"|"+str(lens[i])+"\t"+str(label)+"\n")
cout1.close()

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


#
# deal with R - R
#






'''
cin1= open(com,"r")
index = 0
try:
	while 1:
		line = cin.next().strip()
		line = cin.next().strip()
		array = line.split()
		line = cin.next().strip()
		array1 = line.split()
		t = []
		g = []
		for i in array:
			t.append(i)
		for i in array1:
			g.append(i)
		if len(contigs[array[0]]) > k1:
			temp = contigs[array[0]
			if temp[:k1] in kmerin:
				del kmerin[temp[:k1]]
			else:
				 kmerin[temp[:k1]] = ("U",index,array[0])
		if len(contigs[array[len(array)-1]]) > k1:
			if temp[-k1:] in kmerout:
				del kmerout[temp[-k1:]]
			else:
				kmerout[temp[-k1:]] = ("U",index,array[0])
		else:
			
		index += 1
except:x
	info = sys.exc_info()
	print info[0],info[1]
	cin.close()

cin = open(scf,"r")
try:
	while 1:
		line = cin.next().strip()
		line = cin.next().strip()
		if len(line) >= k1:
			lines.append(line)
			if line[:k1-1] in kmerin:
				del kmerin[line[:k1-1]]
			else:
				 kmerin[line[:k1-1]] = index
			if line[-k1+1:] in kmerout:
				del kmerout[line[-k1+1:]]
			else:
				kmerout[line[-k1+1:]] = index
			index += 1
except:
	info = sys.exc_info()
	print info[0],info[1]
	cin.close()

print "Number of inputed scaffolds and repeat-contigs = " + str(index)
print "Number of in (K-1)mer = " + str(len(kmerin))
print "Number of out (K-1)mer = " + str(len(kmerout))

i = 0
links = [-1 for i in range(len(lines))]

count = 0
index = 0
for out in kmerout:
	if out in kmerin and kmerout[out] != kmerin[out]:
		links[int(kmerout[out])] = kmerin[out]
		count += 1

print "links = " + str(count)
temp = {}
for i in lines:
	temp[index] = i
	index += 1
for i in range(len(lines)):
	line = ''
	flag = links[i]
	if i not in temp:
		continue
	while  flag > -1:
		temp[i] += temp[flag][-k1+1:]
		t = links[flag]
		del temp[flag]
		flag = t
	links[i] = -1
lens = []
maxL = 0
N50 = 0
N90 = 0

sumL = 0
for i in temp:
	lens.append(len(temp[i]))
lens = sorted(lens,reverse=True)
cout = open(prefix+"/scaffold","w")
for i in temp:
	cout.write(">seq_"+str(index)+"\t"+str(len(temp[i]))+"\n")
	cout.write(temp[i]+"\n")
	maxL += len(temp[i]) 
cout.close()
cout = open(prefix+"/lens","w")
for i in lens:
	cout.write(str(i)+"\n")
	sumL += i
	if N50 == 0 and sumL >= maxL/2:
		N50 = i
	if N90 == 0 and sumL >= maxL*0.9:
		N90 = i
cout.close()

print "Longest = " + str(lens[0])
print "N50 = " + str(N50)
print "N90 = " + str(N90)


while k1 <= 51:
	k1 += 2
	lines = []
	kmerin = {}
	kmerout = {}
	index = 0
	cin = open(prefix+"/scaffold","r")
	try:
		while 1:
			line = cin.next().strip()
			line = cin.next().strip()
			if len(line) >= k1:
				lines.append(line)
				if line[:k1-1] in kmerin:
					del kmerin[line[:k1-1]]
				else:
					 kmerin[line[:k1-1]] = index
				if line[-k1+1:] in kmerout:
					del kmerout[line[-k1+1:]]
				else:
					kmerout[line[-k1+1:]] = index
				index += 1
	except:
		info = sys.exc_info()
		print info[0],info[1]
		cin.close()

	print "Number of inputed scaffolds = " + str(index)
	print "Number of in (K-1)mer = " + str(len(kmerin))
	print "Number of out (K-1)mer = " + str(len(kmerout))

	i = 0
	links = [-1 for i in range(len(lines))]

	count = 0
	index = 0
	for out in kmerout:
		if out in kmerin and kmerout[out] != kmerin[out]:
			links[int(kmerout[out])] = kmerin[out]
			count += 1

	print "links = " + str(count)
	temp = {}
	for i in lines:
		temp[index] = i
		index += 1
	for i in range(len(lines)):
		line = ''
		flag = links[i]
		if i not in temp:
			continue
		while  flag > -1:
			temp[i] += temp[flag][-k1+1:]
			t = links[flag]
			del temp[flag]
			flag = t
		links[i] = -1
	lens = []
	maxL = 0
	N50 = 0
	N90 = 0

	sumL = 0
	for i in temp:
		lens.append(len(temp[i]))
	lens = sorted(lens,reverse=True)
	cout = open(prefix+"/scaffold","w")
	for i in temp:
		cout.write(">seq_"+str(index)+"\t"+str(len(temp[i]))+"\n")
		cout.write(temp[i]+"\n")
		maxL += len(temp[i]) 
	cout.close()
	cout = open(prefix+"/lens","w")
	for i in lens:
		cout.write(str(i)+"\n")
		sumL += i
		if N50 == 0 and sumL >= maxL/2:
			N50 = i
		if N90 == 0 and sumL >= maxL*0.9:
			N90 = i
	cout.close()

	print "Longest = " + str(lens[0])
	print "N50 = " + str(N50)
	print "N90 = " + str(N90)
'''

