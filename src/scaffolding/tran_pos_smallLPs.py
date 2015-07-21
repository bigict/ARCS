#!/usr/bin/python

import sys
import os

if len(sys.argv) != 3:
	print "Usage: tran_pos.py workspace iteration"
	os._exit(0)

workspace = sys.argv[1]
if workspace[len(workspace) - 1] != '/':
	workspace += '/'

fin = open(workspace + "scaffold_parameter_" + sys.argv[2], 'r')
edge_num = 0

try:
	while True:
		line = fin.next().strip()
		a, b = line.split()
		if a == 'EDGE_CLUSTER_NUM:':
			edge_num = int(b)
			break
except:
	print "get edge num"
#print "edge num " + str(edge_num)

pos_array = [0 for i in range(edge_num)]
fin.close()

array = []
dir=workspace + "/smallLPResults/"
for name in os.listdir(dir):
	filename = os.path.join(dir, name);
#	print filename
	fin = open(filename, 'r')
	
	fin.next()
	
	line = fin.next().strip()
	(a, b) = line.split()
	#print 'row_num = ' + b
	row_num = int(b)
	line = fin.next().strip()
	(a, b) = line.split()
	#print 'column_num = ' + b
	column_num = int(b)
	
	
	index = 0
	pos = 0
	
	try:
		while True:
			line = fin.next().strip()
			array = line.split()
			if len(array) > 1 and len(array[1]) >2 and array[1][0:2] == "x_" :
				pos = int(array[3])
				index = int(array[1][2:])
				pos_array[index] = pos
	except:
                pass
#		print "read pos file end"
	
	fin.close()
	
fout = open(workspace + "edge_cluster_pos_" + sys.argv[2] , 'w')
	
for i in range(edge_num):
	fout.write(str(pos_array[i]) + "\n")

fout.close()



