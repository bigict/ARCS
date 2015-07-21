#!/usr/bin/python

import sys
import os
import getopt

if len(sys.argv) != 4:
	print "Usage: generateCom.py -i <iter> -D <workspace> -F <scaf_seq>"
	os._exit(0)

print '''[TEST] Generating i+1 component and cdbg using scaf_seq 
		and backup pre component and cdbg\n'''

opts, var = getopt.getopt(sys.argv[1:], "i:D:F:")
for op,value in opts:
	if op == '-i':
		ite = int (value)
	elif op == '-D':
		prefix = value
	elif op == '-F':
		scaf = value

cdbg = prefix + '/cdbg_copy_number.fa'
com  = prefix + '/component_' + str(i+1)

mvcmd = 'mv ' + cdbg + ' ' + cdbg + '_' + str(i)
os.system(mvcmd)
mvcmd = 'mv ' + com + ' ' + com + '_' + str(i)
os.system(mvcmd)

index = 0
fin = open(scaf, 'r')
fout0 = open(cdbg, 'w')
fout1 = open(com , 'w')
try:
    while 1:
        line1 = fin.next().strip()
        line2 = fin.next().strip()
        fout0.write(line1 + "\t1\n" + str(line2) + "\n")
        fout1.write(">component " + str(index) + "\n" + str(index) + "\n\n")
        index += 1
except:
    info = sys.exc_info()
    print info[0],info[1]
fin.close()
fout0.close()
fout1.close()

print "[Done] "
