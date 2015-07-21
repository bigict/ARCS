#!/usr/bin/python

import sys
import os

if len(sys.argv) != 5:
	print "Usage: copy_number_transform.py [input_glpk_output_file] [input_contig_fasta_file] [output_fasta_file_name] [output_component_file]"
	os._exit(0)


print "======================================================================================\n";
print "Step 2. Copy Number\n";
print "======================================================================================\n";
fin = open(sys.argv[1], 'r')
fin.next()

line = fin.next().strip()
(a, b) = line.split()
print 'row_num = ' + b
row_num = int(b)
line = fin.next().strip()
(a, b) = line.split()
print 'column_num = ' + b
column_num = int(b)

fin.next()
fin.next()
fin.next()
fin.next()
fin.next()
fin.next()

count = 0
while count < row_num:
	fin.next().strip()
	count += 1

fin.next().strip()
fin.next().strip()
fin.next().strip()

count = 0
copy_num = {}
array = []
while count < column_num:
	line = fin.next().strip()
	count += 1
	array = line.split()
	
	if len(array) < 4:
		line += '\t'
		line += fin.next().strip()
		array = line.split()

	(a, b) = array[1][2:-1].split(',')
	c =  array[3]
	a = int(a)
	b = int(b)
	if a == 1 or b == row_num:
		continue
	c = float(c)
	if c == 0:
		continue
	if (a, b) not in copy_num:
		copy_num[(a, b)] = c
	else:
		copy_num[(a, b)] += c

#print "count == "
#print count

fin.close()

fin = open(sys.argv[2], 'r')
line = ''
array = []

fout = open(sys.argv[3], 'w')
fout_cmp = open(sys.argv[4], 'w')
count = 0
count_cmp = 0
try:
	while 1:
		line = fin.next().strip()
		array = line.split()
		a = int(array[0])
		b = int(array[1])
		
		if (a, b) in copy_num:
			c = int(round(float(copy_num[(a, b)]),0))
			if c > 0:
				fout.write(">seq\t")
				fout.write(str(c) + '\n')
				fout.write(array[2] + '\n')
				if c == 1:
					fout_cmp.write(">component " + str(count_cmp) + '\n')
					fout_cmp.write(str(count) + '\n')
					fout_cmp.write(' \n')
					count_cmp += 1

				count += 1
		fin.next()
		
except:
	print 'tranform end'
	#print count

fout.close()
fout_cmp.close()
fin.close()

