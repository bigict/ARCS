#!/usr/bin/python
#wirtten by wangbing
import time
import datetime
import os
import commands
import sys
import getopt
import re

usage="Usage:" + sys.argv[0] + " cdbg_copy_num.fa *_ref.fa contig_num_file Kmer_num"

if len(sys.argv) != 5:
    print usage
    os._exit(0);

cdbg_file = sys.argv[1]
ref_file = sys.argv[2]
contig_num_file = sys.argv[3]
predict_time = 0
K = int(sys.argv[4])

fin = open( cdbg_file , 'r' )
fout = open( "cdbg_file.out", 'w')
array = []
count = -1
num_to_contig = []
try:
    while True:
        line = fin.next().strip()
        array = line.split()
        predict_time = int(array[1])
        contig = fin.next().strip()
        count += 1
        fout.write(line + " " + str(count) + "\n")
        fout.write(contig+"\n")
        num_to_contig.append(str(predict_time) + " " + contig)
except:
#    info = sys.exc_info()
#    print info[0],info[1]
#    print "no this num of contig"
    pass

fin.close();

fin = open( ref_file, 'r' )
ref = ''
fin.next()
try:
    while True:
        line = fin.next().strip()
        ref += line
except:
    print "[Info] ref len =", len(ref)
fin.close()

#print ref
kmerFreq = {}
for i in range(len(ref)-K+1):
    kmer = ref[i:i+K]
    if kmer not in kmerFreq:
        kmerFreq[kmer] = 1
    else:
        kmerFreq[kmer] += 1

ref_r = ''
for i in range(len(ref)):
    if ref[len(ref)-1-i] == "A":
        ref_r += "T"
    elif ref[len(ref)-1-i] == "T":
        ref_r += "A"
    elif ref[len(ref)-1-i] == "G":
        ref_r += "C"
    elif ref[len(ref)-1-i] == "C":
        ref_r += "G"
    else:
        ref_r += "N"
assert len(ref) == len(ref_r)
for i in range(len(ref_r)-K+1):
    kmer = ref_r[i:i+K]
    if kmer not in kmerFreq:
        kmerFreq[kmer] = 1
    else:
        kmerFreq[kmer] += 1
#print kmerFreq

def fun(a):
    global num_to_contig
    global kmerFreq
    array = num_to_contig[a].split(" ")
    print "[Info] contig num:", a
    print "[Info] predict freq: " + array[0]
    print "[Info] contig length:", len(array[1])
    print "[Info] contig: " + array[1]
    print "[Info] kmer freq"
    contig = array[1]
    for i in range(len(contig) - K + 1):
        kmer = contig[i:i+K]
        if kmer not in kmerFreq:
            print "0",
        else:
            print kmerFreq[kmer],
    print "\n"

fin = open(contig_num_file, 'r')
try:
    while True:
        line = fin.next().strip()
        num = int(line)
        if num >= len(num_to_contig):
            print "[Info] contig num:", num
            print "[Info] contig over range"
        else:
            fun(num)
except:
    print "[Info] end contig num file"
'''
def Test(a,b,c):
    print "def function..."
    for idx in a:
        print idx
    for i in range(len(b)):
        print b[i]
    return a,b,c,c,c

a = {}
b = [1,2,3]
c = 0
part = Test(a,b,c)
print type(part)
'''

#./check_repeat.py /home/zheng/result/sta/arcs/25mer/cdbg_copy_number.fa /home/zheng/data/Staphylococcus_aureus/Data/original/genome.fasta nums 25

