#!/usr/bin/python

import time
import datetime
import os
import commands
import sys
import getopt
import re

print ""
print "======================================================================================"
print ""
print "                                 ARCS (Q Version 0.9)                               "
print "	    Copyright(c) 2014, Renyu Zhang, Qing Xu, Dongbo Bu. All Rights Reserved.       " 
print ""
print "======================================================================================"
print ""
print time.strftime('%Y-%m-%d',time.localtime(time.time()))
print ""

result = commands.getstatusoutput("glpsol")

if result[0] == 32512:
	print "glpk has not been installed"
	print "please install glpk and run again"
	print "for ubuntu, just sudo apt-get install glpk"
	os._exit(0)

USAGE = "USAGE: ARCS.py [options]\n" + \
 "DESCRIPTION:\n" + \
 "OPTIONS:\n" + \
 "\t-s --configure_file       configure file\n" + \
 "\t-d --workspace            workspace\n" + \
 "\t-K --kmer_size            kmer size\n" + \
 "\t-e --edge_length_cutoff   edge length cutoff\n" + \
 "\t-O --max_overlap          max overlap to detect conflict\n" + \
 "\t-h --help                 help information\n" + \
 "\t-v --version              software version" + \
 "\t-p --CPU                  CPU number to be used"

if len(sys.argv) == 1:
	print USAGE
	os._exit(0)

try:
	longopts = ["help", "version", "configure_file=", "workspace=", "kmer_size=", "edge_length_cutoff=", "max_overlap", "CPU"]
except:
	print "parameter error"
opts, var = getopt.getopt(sys.argv[1:], "s:d:K:e:hO:vp:", longopts)

configure_file = ''
workspace = ''
kmer_size = 27
edge_length_cutoff = kmer_size
max_overlap = 200
cpu_num = 8

if len(opts) == 0:
	print USAGE
	os._exit(0)
for pair in opts:
	if pair[0] == '-h' or pair[0] == '--help':
		print USAGE
		os._exit(0)
	elif pair[0] == '-s' or pair[0] == '--configure_file':
		configure_file = pair[1]
	elif pair[0] == '-d' or pair[0] == '--workspace':
		workspace = pair[1]
	elif pair[0] == '-K' or pair[0] == '--kmer_size':
		kmer_size = int(pair[1])
		if kmer_size < 1 or kmer_size > 96:
			print "kmer size must be in [1, 96]"
			print "system exit"
			os._exit(0)
	elif pair[0] == '-e' or pair[0] == '--edge_length_cutoff':
		edge_length_cutoff = int(pair[1])
	elif pair[0] == '-p' or pair[0] == '--CPU':
		cpu_num = int(pair[1])
	elif pair[0] == '-O' or pair[0] == '--max_overlap':
		max_overlap = int(pair[1])
	elif pair[0] == '-v' or pair[0] == '--version':
		print "Version 0.9: released on July 2th, 2014"
		os._exit(0)
	else:
		print "no option -" + pair[0]
		print USAGE
		os._exit(0)

if workspace == '':
	print "please specify a workspace"
	os._exit(0)

if configure_file == '':
	print "please specify a configure file"
	os._exit(0)

word = ''
if not os.path.exists(workspace):
	print "directory " + workspace + " does not exist"
	print "creat this directory yes/no"
	word = sys.stdin.readline().strip()
	word = word.lower()
	if word == 'y' or word == 'yes':
		os.mkdir(workspace)
	else:
		print "system exit"
		os._exit(0)

if workspace[len(workspace) - 1] == '/':
	workspace = workspace[0 : len(workspace) - 1]
workspace = os.path.abspath(workspace)
print "path : " + workspace

start = datetime.datetime.now()

fin = open(configure_file, 'r')
array = []
cha = ' '
lib_list = []
insert_size = []
link_quality_percent = []
pair_kmer_cutoff = []
pair_reads_cutoff = []

q1 = ''
q2 = ''

try:
	while True:
		line = fin.next().strip()
		for i in range(len(line)):
			if line[i] != ' ':
				cha = line[i]
		if cha == '#':
			continue
		line = line.strip()
		array = re.split('[:=]', line)
		if array[0] == 'q1':
			q1 = array[1]
			q1 = q1.strip()
		if array[0] == 'q2':
			q2 = array[1]
			q2 = q2.strip()
			lib_list.append((q1, q2))
		if array[0] == 'INSERT_SIZE':
			insert_size.append(array[1].strip())
		if array[0] == 'LINK_QUALITY_PERCENT':
			link_quality_percent.append(array[1].strip())
		if array[0] == 'PAIR_KMER_CUTOFF':
			pair_kmer_cutoff.append(array[1].strip())
		if array[0] == 'PAIR_READS_CUTOFF':
			pair_reads_cutoff.append(array[1].strip())
		

except:
	print "pass configure file"
'''
if kmer_size > 0 and kmer_size <= 32:
	contiging_cmd = 'contiging/32-mer/contiging -s ' + configure_file + ' -d ' + workspace + ' -K ' + str(kmer_size)
elif kmer_size > 32 and kmer_size <= 64:	
	contiging_cmd = 'contiging/64-mer/contiging -s ' + configure_file + ' -d ' + workspace + ' -K ' + str(kmer_size)
elif kmer_size <= 96:	
	contiging_cmd = 'contiging/96-mer/contiging -s ' + configure_file + ' -d ' + workspace + ' -K ' + str(kmer_size)

if os.system(contiging_cmd) != 0:
	os._exit(1)

glpsol_cmd = 'glpsol --mincost ' + workspace + '/min_cost_flow.DIMACS -o ' + workspace + '/min_cost.out  > /dev/null'
if os.system(glpsol_cmd) != 0:
	os._exit(1)

copy_num_cmd = 'contiging/copy_number_transform.py ' + workspace + '/min_cost.out ' + workspace + '/condensed_de_bruijn_graph_after_trimming.data ' + workspace + '/cdbg_copy_number.fa ' + workspace + '/component_0' 

if os.system(copy_num_cmd) != 0:
	os._exit(1)


print "link quality size " + str(len(link_quality_percent))
'''

for i in range(len(lib_list)):
	print "............ iter " + str(i + 1)
	ele = lib_list[i]

	scaffolding_cmd = 'Rscaffolding/Rscaffolding -d ' + workspace + ' -K ' + str(kmer_size) + ' -c ' + workspace + '/cdbg_copy_number.fa -e ' + str(edge_length_cutoff) + ' -1 ' + ele[0] + ' -2 ' + ele[1] + ' -L ' + insert_size[i] + ' -P ' + link_quality_percent[i] + ' -p ' + str(cpu_num) + ' -i ' + str(i) + ' -r ' + pair_kmer_cutoff[i] + ' -R ' + pair_reads_cutoff[i] 
	if os.system(scaffolding_cmd) != 0:
        	os._exit(1)

	glpsol_cmd = 'glpsol --math ' + workspace + '/position_lp_' + str(i) + '.math -o ' + workspace + '/position_lp_' + str(i) + '.out  > /dev/null'
	if os.system(glpsol_cmd) != 0:
		os._exit(1)
	tran_pos_cmd = 'scaffolding/tran_pos.py ' + workspace + ' ' + str(i)
	if os.system(tran_pos_cmd) != 0:
		os._exit(1)
        
	remove_repeats_cmd = 'remove_repeats/remove_repeats -d ' + workspace + ' -O ' + str(max_overlap) + ' -K ' + str(kmer_size) + ' -i ' + str(i)
	if os.system(remove_repeats_cmd) != 0:
		os._exit(1)

	gap_filling_cmd = 'gap_filling/gap_filling -s scaffold_parameter_'+str(i)+' -K ' + str(kmer_size) + ' -O ' + str(kmer_size - 10) + ' -c cdbg_copy_number.fa -l component_' + str(i+1) + ' -d ' + workspace + ' -i condensed_de_bruijn_graph_before_trimming.data'
	if os.system(gap_filling_cmd) != 0:
		os._exit(1)
	
        if int(i+1) in range(len(lib_list)):
		generate_com = './generateCom.py -i ' + str(i) + ' -F ' + workspace + '/' + str(kmer_size) + 'mer.scaf_seq_with_gaps ' + ' -D ' + workspace + '/'
		if os.system(gap_filling_cmd) != 0:
			os._exit(1)
	
end = datetime.datetime.now()
print 'total running time is ' +  str((end - start).seconds) + ' seconds'
	
	
#i = 0 
#ele = lib_list[i]
#component_cmd = './Rcomponent.py ' + str(i) + ' ' +  workspace
#Scomponent_cmd = './Scomponent.py ' + str(i) + ' ' +  workspace + ' ' + str(kmer_size)+"mer.scaf_seq_with_gaps"
#if os.system(component_cmd) != 0:
#	os._exit(1)



'''
print "link quality size " + str(len(link_quality_percent))


for i in range(len(lib_list)):
	print "............ iter " + str(i + 1)
	ele = lib_list[i]

        scaffolding_cmd = 'scaffolding/scaffolding -d ' + workspace + ' -K ' + str(kmer_size) + ' -c ' + workspace + '/cdbg_copy_number.fa -e ' + str(edge_length_cutoff) + ' -1 ' + ele[0] + ' -2 ' + ele[1] + ' -L ' + insert_size[i] + ' -P ' + link_quality_percent[i] + ' -p ' + str(cpu_num) + ' -i ' + str(i) + ' -r ' + pair_kmer_cutoff[i] 
 	
        if os.system(scaffolding_cmd) != 0:
   	        os._exit(1)

	glpsol_cmd = 'glpsol --math ' + workspace + '/position_lp_' + str(i) + '.math -o ' + workspace + '/position_lp_' + str(i) + '.out  > /dev/null'
	if os.system(glpsol_cmd) != 0:
		os._exit(1)
	tran_pos_cmd = 'scaffolding/tran_pos.py ' + workspace + ' ' + str(i)
	if os.system(tran_pos_cmd) != 0:
		os._exit(1)
        
	remove_repeats_cmd = 'remove_repeats/remove_repeats -d ' + workspace + ' -O ' + str(max_overlap) + ' -K ' + str(kmer_size) + ' -i ' + str(i)
	if os.system(remove_repeats_cmd) != 0:
		os._exit(1)




gap_filling_cmd = 'gap_filling/gap_filling -s scaffold_parameter_'+str(len(insert_size)-1)+' -K ' + str(kmer_size) + ' -O ' + str(kmer_size - 10) + ' -c cdbg_copy_number.fa -l component_' + str(len(insert_size)) + ' -d ' + workspace + ' -i condensed_de_bruijn_graph_before_trimming.data'

if os.system(gap_filling_cmd) != 0:
	os._exit(1)
'''
