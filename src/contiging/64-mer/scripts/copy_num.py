#!/usr/bin/python

fin = open('min_out.txt', 'r')

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

print "count == "
print count
fin.close()

file_name = 'condensed_de_bruijn_graph.txt'

fin = open(file_name, 'r')
line = ''
array = []
kmer = ''

fout = open('cdbg_cp_num.txt', 'w')
count = 0
try:
	while 1:
		line = fin.next().strip()
		array = line.split()
		a = int(array[0])
		b = int(array[1])
		
		if (a, b) in copy_num:
			count += 1
			c = int(round(float(copy_num[(a, b)]),0))
			if c > 0:
				fout.write(">seq\t")
				fout.write(str(c) + '\n')
				fout.write(array[2] + '\n')

		fin.next()
		
except:
	print 'end !!'
	print count

fout.close()
fin.close()

