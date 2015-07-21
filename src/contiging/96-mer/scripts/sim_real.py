#!/usr/bin/python

K = 25

ref_file = 'data/ecoli_ref.fa'

ref_in = open(ref_file, "r")
ref = ''
ref_in.next()

def reversecomplement( read ):
	rcread = ''
	#print read
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

try:
	while 1:
		ref += ref_in.next().strip()
except:
	print 'read ref end'
ref_in.close()

ref_r = reversecomplement(ref)

kmerFreq = {}

for i in range(len(ref)-K+1):
	kmer = ref[i:i+K]
	if kmer in kmerFreq:
		kmerFreq[kmer] += 1
	else:
		kmerFreq[kmer] = 1
for i in range(len(ref_r)-K+1):
	kmer = ref_r[i:i+K]
	if kmer in kmerFreq:
		kmerFreq[kmer] += 1
	else:
		kmerFreq[kmer] = 1

file_name = 'condensed_de_bruijn_graph.txt'

fin = open(file_name, 'r')
line = ''
array = []
kmer = ''

fout = open('cdbg_cp_num.txt', 'w')

i = 0

try:
	while 1:
		line = fin.next().strip()
		array = line.split()
		line = array[2]
		kmer = line[0:K]
		
		fin.next()
		
		if kmer in kmerFreq:
			fout.write(">seq" + str(i) + "\t" + str(kmerFreq[kmer]) + "\n")
		else:	
			fout.write(">seq" + str(i) + "\t0\n")

		fout.write(line + '\n')

		i += 1
			
except:
	print 'end !!'

fout.close()
fin.close()

