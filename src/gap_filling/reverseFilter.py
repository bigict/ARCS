#!/usr/bin/python
Usage = '''
	$python reverseFilter.py <workspace> <cdbg> <component> 
	'''

import sys
import os

if len(sys.argv) < 4:
	print Usage
	os._exit(0)

path = sys.argv[1]
cdbgFile = sys.argv[2]
comFile = sys.argv[3]

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


count = 0
reverse = {}	# i to j

cdbg = {}
cdbgs = {}
lens = {}
fin = open(cdbgFile, 'r')
count = 0;
try:
	while 1:
		line = fin.next()
		array = line.split()
		if int(array[1]) > 1:
			count += 1
			line = fin.next().strip();
			continue
		# rep[count] = int(array[1])
		line = fin.next().strip();
		# array = line.split()
		cdbg[line] = count
		cdbgs[count] = line
		lens[count] = len(line)
		count += 1
except:
	info = sys.exc_info()
	print info[0],info[1]
	fin.close()

for ctg in cdbg:
	pos = cdbg[ctg]
	if lens[pos] < 100:
		continue
	if pos in reverse:
		continue
	re = reversecomplement(ctg)
	if re in cdbg:
		reverse[pos] = cdbg[re]
		reverse[cdbg[re]] = pos

print "[info] num of unique ctg = ",len(cdbg)
print "[info] num of paired ctg = ",len(reverse)
fout = open(path+"/misspaired.fa",'w')
for ctg in cdbg:
	pos = cdbg[ctg]
	if pos not in reverse:
		fout.write(">ctg\n"+ctg+"\n"+reversecomplement(ctg)+"\n")
fout.close()

count = 0
scflens = {}
sumL = 0
tmp = []
scfs = []	# [[ids]] ids of ctgs
ctg2scf = {}	#id of ctg 2 id of scf
original = []

fin = open(comFile,'r')
try:
	while 1:
		line1 = fin.next()
		line1 = fin.next().strip()
		ids = line1.split()
		line2 = fin.next().strip()
		gaps = line2.split()
		original.append((line1, line2))
		scfs.append(ids)
		l = 0
		for i in range(len(ids)):
			ctg2scf[int(ids[i])] = count
			l += lens[int(ids[i])]
		for i in gaps:
			l += int(i)
		sumL += l
		scflens[count] = l
		tmp.append((count,l))
		count += 1
except:
	info = sys.exc_info()
	print info[0],info[1]
	fin.close()
print "[info] num of scfs = ",len(scfs)

sortedLens = [x[0] for x in sorted(tmp,key=lambda x:x[1], reverse=True)]
fout = open(path+"/scflens_predict_0.data",'w')
for l in sortedLens:
    fout.write(str(scflens[l])+"\n")
fout.close()

#print scflens[sortedLens[0]],scflens[sortedLens[1]]
tmpl = 0
print "[info] sumG = ",sumL
for index in sortedLens:
	tmpl += scflens[index]
	if tmpl	>= sumL/2:
		print "[info] N50 = ",scflens[index]
		break

del tmp
flag = [ {} for i in scfs]
for index in sortedLens:
	for ctg in scfs[index]:
		ctg = int(ctg)
		if ctg in reverse:
			re = reverse[ctg]
		else:
			continue
		if re in ctg2scf:
			sid = ctg2scf[re]
			if sid in flag[index]:
				flag[index][sid] += 1
			else:
				flag[index][sid] = 1

num = 0
sumL = 0
newscfs = []
dup = set()
dflag = [0 for i in scfs]
for index in sortedLens:
	if dflag[index] == 1:
		continue
	newscfs.append(index)
	sumL += scflens[index]
	for sid in flag[index]:
		if flag[index][sid] <= 1:
			continue
		if dflag[sid] == 1:
			num += 1
			dup.add(sid)
		dflag[sid] = 1
print "[info] duplication = ",num
print "[info] remain scfs = ",len(newscfs)

#fout = open("dump.lens",'w')
#for index in dup:
#	fout.write(str(scflens[index])+"\n")
#fout.close()

count = 0
fout = open(path+"/component_last",'w')
for index in newscfs:
	fout.write(">component\t"+str(count)+"\n")
	fout.write(original[index][0]+"\n"+original[index][1]+"\n")
	count += 1
fout.close()
fout = open(path+"/scflens_predict_1.data",'w')
for index in newscfs:
    fout.write(str(scflens[index])+"\n")
fout.close()

tmpl = 0
print "[info] sumG = ",sumL
for index in newscfs:
	tmpl += scflens[index]
	if tmpl	>= sumL/2:
		print "[info] N50 = ",scflens[index]
		break

