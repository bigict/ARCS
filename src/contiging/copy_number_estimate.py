#!/usr/bin/python

import sys
import os

if len(sys.argv) != 5:
	print "Usage: copy_number_estimate.py [contig file] [parameter file] [output_fasta_file_name] [component0]"
	os._exit(0)

fin = open(sys.argv[2], 'r')
array = []

lam = 129.0
K = 31

try:
	while True:
		line = fin.next().strip()
		array = line.split("=")
		#print array
		if array[0] == 'lambda':
			lam = float(array[1])
		if array[0] == 'K':
			K = int(array[1])
                        print "K = ",str(K)

except:
	#print "read parameter end"

        fin.close()


contig = ''
cov = []
sum = 0
index = 0;
num = 0
count = 0

fin = open(sys.argv[1], 'r')

fout = open(sys.argv[3], 'w')
fout1 = open(sys.argv[4], 'w')

try:
	while True:
		line = fin.next().strip()
		array = line.split()
		contig = array[2]
		line = fin.next().strip()
		cov = line.split()
		if len(cov) != len(contig) - K + 1:
                        print "[info] \t"+contig,count
			print "length not ok"
			count += 1
                        continue
                        #sys.exit(0)
		sum = 0
		for i in range(len(cov)):
			sum += int(cov[i])
		copy_num = int(round((sum / len(cov))/lam))
		if copy_num >= 0:
			fout.write(">seq_" + str(index) + " \t" + str(copy_num) + '\n')
			fout.write(contig + '\n')
                        if copy_num <= 1:
			        fout1.write(">component\t" + str(num) + '\n')
			        fout1.write(str(index) + '\n\n')
                                num += 1
                        
                        index += 1
                count += 1

except:
        fin.close()
        fout.close()
        fout1.close()
	print "\tNumber of condensed edges = ",str(count)
	print "\tNumber of Contigs = ",str(index)
	print "\tNumber of Uniques = ",str(num)

