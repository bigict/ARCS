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
	fin.next()
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
fout = open('predict.txt', 'w')

count = 0
for pair in copy_num:
	fout.write(str(pair[0]) + '\t' + str(pair[1]) + '\t' + str(copy_num[pair]) + '\n')
	count += 1

print count
fout.close()
