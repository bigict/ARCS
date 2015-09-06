#!/usr/bin/python

import time
from datetime import datetime
import os
import commands
import sys
import getopt
import re

print ""
print "======================================================================================"
print ""
print "                                 ARCS (Q Version 0.9)                               "
print "        Copyright(c) 2014, Renyu Zhang, Qing Xu, Dongbo Bu. All Rights Reserved.       " 
print ""
print "======================================================================================"
print ""
print datetime.now()
print ""

result = commands.getstatusoutput("glpsol --help")

if result[0] != 0:
    print "glpk has not been installed"
    print "please install glpk and run again"
    print "for ubuntu, just sudo apt-get install glpk"
    sys.exit(1)

USAGE = """"USAGE: ARCS.py [options]
 DESCRIPTION:
 OPTIONS:
    -s --configure_file       configure file
    -d --workspace            workspace
    -K --kmer_size            kmer size
    -e --edge_length_cutoff   edge length cutoff
    -O --max_overlap          max overlap to detect conflict
    -h --help                 help information
    -v --version              software version
    -E --kmer_filter          filter low quality kmers
    -p --CPU                  CPU number to be used
"""

if len(sys.argv) == 1:
    print USAGE
    os._exit(0)

try:
    longopts = ["help", "version", "configure_file=", "workspace=", "kmer_size=", "edge_length_cutoff=", "max_overlap", "CPU","kmer_filter"]
except:
    print "parameter error"
opts, var = getopt.getopt(sys.argv[1:], "s:d:K:e:hO:vp:E:c:", longopts)

configure_file = ''
workspace = ''
kmer_size = 27
edge_length_cutoff = kmer_size
max_overlap = 200
cpu_num = 8
kmer_filter = 0

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
    elif pair[0] == '-E' or pair[0] == '--kmer_filter':
        kmer_filter = int(pair[1])
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

start = datetime.now()

fin = open(configure_file, 'r')
array = []
cha = ' '
lib_list = []
insert_size = []
link_quality_percent = []
pair_kmer_cutoff = []
pair_reads_cutoff = []
contig_lengths_cutoff = []

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
        if array[0] == 'EDGE_LENGTH_CUTOFF':
            contig_lengths_cutoff.append(array[1].strip())
        

except:
    print "pass configure file"

Path = sys.path[0] + "/"
print Path

if kmer_size > 0 and kmer_size <= 32:
    if kmer_filter:
        contiging_cmd = Path + 'contiging/contiging -s ' + configure_file + ' -d ' + workspace + ' -K ' + str(kmer_size) + ' -E ' + ' -p ' + str(cpu_num)
    else:
        contiging_cmd = Path + 'contiging/contiging -s ' + configure_file + ' -d ' + workspace + ' -K ' + str(kmer_size) + ' -p ' + str(cpu_num)
elif kmer_size > 32 and kmer_size <= 64:    
    contiging_cmd = Path + 'contiging/contiging -s ' + configure_file + ' -d ' + workspace + ' -K ' + str(kmer_size)
elif kmer_size <= 96:    
    contiging_cmd = Path + 'contiging/contiging -s ' + configure_file + ' -d ' + workspace + ' -K ' + str(kmer_size)

print '-------------------------------------'
print contiging_cmd
print '-------------------------------------'

if os.system(contiging_cmd) != 0:
    os._exit(1)

copy_num_cmd = 'cat ' + workspace + '/condensed_de_bruijn_graph_after_trimming.data | ' + Path + 'contiging/copy_number_estimate.py -p' + workspace + '/contig_parameter ' + workspace + '/cdbg_copy_number.fa ' + workspace + '/component_0' 

print '-------------------------------------'
print copy_num_cmd
print '-------------------------------------'

if os.system(copy_num_cmd) != 0:
    os._exit(1)


print "link quality size " + str(len(link_quality_percent))


for i in range(len(lib_list)):
    print "............ iter " + str(i + 1)
    ele = lib_list[i]

    if len(contig_lengths_cutoff) <= i:
        contig_lengths_cutoff.append(edge_length_cutoff)
    #scaffolding_cmd = Path + 'scaffolding/scaffolding -d ' + workspace + ' -K ' + str(kmer_size) + ' -c ' + workspace + '/cdbg_copy_number.fa -e ' + str(contig_lengths_cutoff[i]) + ' -1 ' + ele[0] + ' -2 ' + ele[1] + ' -L ' + insert_size[i] + ' -P ' + link_quality_percent[i] + ' -p ' + str(cpu_num) + ' -i ' + str(i) + ' -r ' + pair_kmer_cutoff[i] + ' -R ' + pair_reads_cutoff[i] 
    scaffolding_cmd = Path + 'scaffolding/scaffolding -d ' + workspace + ' -K ' + str(kmer_size) + ' -C ' + workspace + '/cdbg_copy_number.fa -f ' + workspace + '/component_' + str(i)  + ' -e ' + str(contig_lengths_cutoff[i]) + ' -1 ' + ele[0] + ' -2 ' + ele[1] + ' -L ' + insert_size[i] + ' -P ' + link_quality_percent[i] + ' -p ' + str(cpu_num) + ' -i ' + str(i) + ' -r ' + pair_kmer_cutoff[i] + ' -R ' + pair_reads_cutoff[i] 
    print '-------------------------------------'
    print scaffolding_cmd
    print '-------------------------------------'
     
    if os.system(scaffolding_cmd) != 0:
        os._exit(1)

    #change lp to smallLPs by wangbing   
    glpsol_cmd = Path + 'divideLP/runLP.sh ' +  workspace + '/position_lp_' + str(i) + '.math ' + workspace + '/smallLPs/ ' + workspace + '/smallLPResults/'
    print '-------------------------------------'
    print glpsol_cmd
    print '-------------------------------------'
    if os.system(glpsol_cmd) != 0:
        os._exit(1)

    tran_pos_cmd = Path + 'scaffolding/parse_glpk_results.py -d ' + workspace + ' -i ' + str(i) + ' -p ' + os.path.join(workspace, 'scaffold_parameter_%d' % (i))
    print '-------------------------------------'
    print tran_pos_cmd
    print '-------------------------------------'
    if os.system(tran_pos_cmd) != 0:
        os._exit(1)
    #end change
    remove_repeats_cmd = Path + 'remove_repeats/remove_repeats -d ' + workspace + ' -O ' + str(max_overlap) + ' -K ' + str(kmer_size) + ' -i ' + str(i)
    print '-------------------------------------'
    print remove_repeats_cmd 
    print '-------------------------------------'

    if os.system(remove_repeats_cmd) != 0:
        os._exit(1)
        #kmer_size -= 2

reverse_filter_cmd = Path + 'gap_filling/reverse_filter.py ' + workspace + ' ' + workspace + '/cdbg_copy_number.fa ' + workspace + '/component_' + str(len(insert_size))
print '-------------------------------------'
print reverse_filter_cmd
print '-------------------------------------'
if os.system(reverse_filter_cmd) != 0:
	os._exit(1)

gap_filling_cmd = Path + 'gap_filling/gap_filling -s scaffold_parameter_0 -K ' + str(kmer_size) + ' -O ' + str(kmer_size - 10) + ' -c cdbg_copy_number.fa -l component_last -d ' + workspace + ' -i condensed_de_bruijn_graph_before_trimming.data'
print '-------------------------------------'
print gap_filling_cmd
print '-------------------------------------'

end = datetime.now()
print 'total running time is ' +  str((end - start).seconds) + ' seconds'

