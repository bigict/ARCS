#!/usr/bin/python
Usage='''
Author: Zheng QuanGang
Function: Copy number accuracy using kmer-cov and de Bruijn graph

$python copy_number_accuracy.py -P <contig parameter file> -D <workspace> -B <condensed_de_bruijn_graph_after_trimming> -C <cdbg_copy_number.fa>
'''
import sys,os
import getopt

opts, args = getopt.getopt(sys.argv[1:],"P:B:C:D:")

if len(opts) != 4:
    print Usage
    sys.exit()

for op, value in opts:
    if op == "-P":
        para = value
    elif op == "-B":
        debruijn = value
    elif op == "-C":
        cdbg = value
    elif op == "-D":
        outPrefix = value
    else:
        print "Parameter error!"

fin = open(para, 'r')
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
except:
	print "K = " + str(K)
	print "lambda = " + str(lam)
	print "read parameter end"
fin.close()

if os.path.exists(outPrefix+"/copy_number")==False:
    os.system("mkdir "+outPrefix+"/copy_number")
outPrefix = outPrefix+"/copy_number/"

fin = open(cdbg, 'r')
line = ''
array = []
copy_num = 0
count = 0
uniq = 0
multi = 0

try:
	while 1:
		line = fin.next().strip()
		array = line.split()
		copy_num = int(array[1])
		line = fin.next().strip()
		if copy_num == 1:
			uniq += 1
		else:
			multi += 1
		count += 1			
except:
	fin.close()
        info = sys.exc_info()
        print info[0],info[1]
print "\n[TEST]\tCopy numbers of contigs In original CDBG file:"
print '\tContigs in all = ' + str(count) 
print '\tUnique contigs = ' + str(uniq)
print '\tRepeat contigs = ' + str(multi)


out = open(outPrefix+"/cdbd_copy_number","w")
fout = open(outPrefix+'/component', 'w')
fin = open(debruijn,"r")
innode = {}
outnode = {}
count = 0
uniq = 0
multi = 0
index = 0
cindex = 0
try:
	while 1:
		line =  fin.next().strip()
		array = line.split()
		if array[0] in outnode:
			outnode[array[0]] += 1
		else:
			outnode[array[0]] = 1
			innode[array[0]] = 0
		if array[1] in innode:
			innode[array[1]] += 1
		else:
			innode[array[1]] = 1
			outnode[array[1]] = 0
		line =  fin.next().strip()
except:
	fin.close()
        info = sys.exc_info()
        print info[0],info[1]
fin = open(debruijn,"r")	
try:
	while True:
		count += 1
		copy_num = 0
		line = fin.next().strip()
		array = line.split()
		contig = array[2]
                #if len(contig) < K + 2:
		#	line = fin.next().strip()
			#print "continue"
		#	continue
		line = fin.next().strip()
		cov = line.split()
		if len(cov) != len(contig) - K + 1:
			#print len(cov),len(contig),K
			print "\tsomething wrong happened : len(cov) != len(contig) - K + 1"
			sys.exit(0)
		
		if array[0] in innode and array[1] in outnode:
                        if int(innode[array[0]]) > 1 or int(outnode[array[1]]) > 1 or len(contig) < 4*K:
				sumc = 0
				for i in range(len(cov)):
					sumc += int(cov[i])
				copy_num = int(round((sumc / len(cov))/ (lam*1.15) ))
			else:
				copy_num = 1
		if copy_num > 0:
			out.write(">seq_" + str(index) + "_" + str(len(contig)) + "\t" + str(copy_num) + '\n')
			out.write(contig + '\n')
			if copy_num == 1:
				uniq += 1
				fout.write(">component " + str(cindex) + "\n" + str(index) + "\n\n")
				index += 1
				cindex += 1
			else:
				index += 1
				multi += 1
except:
	fin.close()
	out.close()
	fout.close()
        info = sys.exc_info()
        print info[0],info[1]
        print "[DONE]"
print "\n[TEST]\tCopy numbers of contigs In new CDBG file:"
print '\tContigs in all(>K+2) = ' + str(index) + "/" + str(count) 
print '\tUnique contigs = ' + str(uniq)
print '\tRepeat contigs = ' + str(multi)

