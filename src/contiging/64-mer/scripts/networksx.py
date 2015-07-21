#!/usr/bin/python

import networkx

G = networkx.DiGraph()

filename = 'min_cost_flow.DIMACS'
fin = open(filename, 'r')

line = fin.next().strip()
array = line.split()

line = fin.next().strip()
array = line.split()
G.add_node(array[1], demand = - int(array[2]))

line = fin.next().strip()
array = line.split()
G.add_node(array[1], demand = -int(array[2]))

try:
	while 1:
		line = fin.next().strip()
		array = line.split()
		print array
		G.add_edge(array[1], array[2], weight = float(array[5]), capacity = float(array[4]))
except:
	print "end"	
flowDict = networkx.min_cost_flow(G)
print flowDict
