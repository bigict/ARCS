#! /usr/bin/python

Usage = '''
Author : Zheng QuanGang
Function : 

$python coverage.py -K <kmer> 
		   -D <workspace>
		   -S <seq>
		   <reads 1> <reads 2>
'''

import sys,getopt
import os,math

    
contigs = []
covs = []
name = []
index = 0
kmers = {}
line = ''
mer = ''
seq_file = ""
K = 31
prefix = "./"

# Function : Revese sequence
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

# Reading seq file
def readseq():
    global seq_file
    global kmers
    global contigs
    global name
    global index
    global K
    tmpdel = []
    
    print seq_file
    cin = open(seq_file,"r")
    try:
    	while 1:
    		line = cin.next().strip()
    		name.append(line)
    		line = cin.next().strip()
    		contigs.append(line)
                cov = []
    		for i in range(len(line)-K+1):
    		    tmp = line[i:i+K]
    		    if tmp in kmers:
    			kmers[tmp].append((index,i))
    		    else:
                        kmers[tmp] = []
    			kmers[tmp].append((index,i))
                    cov.append(0)
                covs.append(cov)
    		index += 1
    except:
    	info = sys.exc_info()
    	print info[0],info[1]
    	cin.close()

    print "Number of kmers = " + str(len(kmers))
    

# Reading reads file and trainning insert-size
def mapreads():
    cin1 = open(rds1,"r")
    cin2 = open(rds2,"r")
    count = 0
    try:
    	while 1:
    		temp1 = cin1.next().strip()
    		temp1 = cin1.next().strip()
    		temp2 = cin2.next().strip()
    		temp2 = cin2.next().strip()
    		R1 = reversecomplement(temp1)
    		R2 = reversecomplement(temp2)
    		minl = min(len(temp1),len(temp2))
    		for i in range(minl-K+1):
    			kmer1 = temp1[i:i+K]
    			kmer2 = temp2[i:i+K]
                        if kmer1 in kmers:
                            for p in kmers[kmer1]:
                                i = p[0]
                                j = p[1]
                                covs[i][j] += 1
                        if kmer2 in kmers:
                            for p in kmers[kmer2]:
                                i = p[0]
                                j = p[1]
                                covs[i][j] += 1
    		for i in range(minl-K+1):
    			kmer1 = R1[i:i+K]
    			kmer2 = R2[i:i+K]
                        if kmer1 in kmers:
                            for p in kmers[kmer1]:
                                i = p[0]
                                j = p[1]
                                covs[i][j] += 1
                        if kmer2 in kmers:
                            for p in kmers[kmer2]:
                                i = p[0]
                                j = p[1]
                                covs[i][j] += 1
    		temp1 = cin1.next().strip()
    		temp1 = cin1.next().strip()
    		temp2 = cin2.next().strip()
    		temp2 = cin2.next().strip()
                print "\r\treads = "+str(count),
                count += 1
    except:
        cin1.close()
        cin2.close()

    fout = open(prefix+"/seqcov.fa","w")
    for s in range(len(covs)):
        cov = covs[s]
        fout.write(name[s]+"\n")
        for i in range(len(cov)):
            fout.write(str(cov[i])+" ")
        fout.write("\n")
    fout.close()
            
    print "[Done] mapping reads file and output sequence coverage"


def main():
    global seq_file
    global prefix
    global rds1
    global rds2
    global K
    opts, args = getopt.getopt(sys.argv[1:],"I:K:C:S:D:M:R:")
    if len(opts) < 3:
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
    	elif op == '-S':
    	    seq_file = value
    	elif op == '-D':
    	    prefix = value
    	else:
    	    print Usage
    	    sys.exit()
    print "seq_file",seq_file

    readseq()
    mapreads()

if __name__ == "__main__":
    main()
