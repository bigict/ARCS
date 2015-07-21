#! /usr/bin/python

Usage = '''
Author : Zheng QuanGang
Function : 

$python pairedkmer.py	<seq.fa>
                        <workspace>
                        <read 1> <read 2>
                        //
'''

import sys

if len(sys.argv) != 5:
	print Usage
	sys.exit()
else:
	prefix = sys.argv[2]
	r1 = sys.argv[3]
	r2 = sys.argv[4]
        seq = sys.argv[1]
	
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
kmers = []
name = []
def readseq():
	global contigs
	global kmers
	contigs = []
        kmers = []
	cin = open(seq,"r")
	try:
		while 1:
			line =cin.next().strip()
                        name.append(line)
			line = cin.next().strip()
			contigs.append(line)
	except:
		cin.close()

def kmerhash():
	global contigs
        global kmers

        for i in range(len(contigs)):
            tmp = {}
            kmer = ""
            for j in range(len(contigs[i])-30):
                kmer = contigs[i][j:j+31]
                tmp[kmer] = []
                tmp[kmer].append(j)
                tmp[kmer].append(0)
            kmers.append(tmp)

def mapcount():
	global r1
        global r2
        global kmers
        global name
        K = 31
        cin1 = open(r1,"r")
        cin2 = open(r2,"r")
        try:
        	while 1:
        		temp1 = cin1.next().strip()
        		temp1 = cin1.next().strip()
        		temp2 = cin2.next().strip()
        		temp2 = cin2.next().strip()
        		R2 = reversecomplement(temp2)
        		minl = min(len(temp1),len(temp2))
                        assert minl == 36
        		for i in range(minl-30):
        			kmer1 = temp1[i:i+K]
        			kmer2 = R2[i:i+K]
                                for j in range(len(kmers)):
                                        if kmer1 in kmers[j] and kmer2 in kmers[j]:
					    kmers[j][kmer1][1] += 1
					    kmers[j][kmer2][1] += 1
        
        		R1 = reversecomplement(temp1)
        		minl = min(len(temp1),len(temp2))
                        assert minl == 36
        		for i in range(minl-30):
        			kmer1 = R1[i:i+K]
        			kmer2 = temp2[i:i+K]
                                for j in range(len(kmers)):
                                        if kmer1 in kmers[j] and kmer2 in kmers[j]:
					    kmers[j][kmer1][1] += 1
					    kmers[j][kmer2][1] += 1
        
        		temp1 = cin1.next().strip()
        		temp1 = cin1.next().strip()
        		temp2 = cin2.next().strip()
        		temp2 = cin2.next().strip()
        except:
            info = sys.exc_info()
            print info[0],info[1]
        cin1.close()
        cin2.close()
        
	fout = open(prefix+"/paired.fa","w")
        for i in range(len(kmers)):
            fout.write(name[i]+"_paired"+"\n")
            tmp = {}
            for j in kmers[i]:
                tmp[kmers[i][j][0]] = kmers[i][j][1]
            assert len(kmers[i]) == len(contigs[i])-30
            for j in range(len(tmp)):
                fout.write(str(tmp[j])+" ")
            fout.write("\n")

readseq()
kmerhash()
mapcount()
