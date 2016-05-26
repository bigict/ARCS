#!/usr/bin/python
#
import sys
import os
from itertools import ifilter,imap
import string
import re

def make_complement(read):
    rcread = ''
    for c in read:
        if c == 'A':
            rcread += 'T'
        elif c == 'C':
            rcread += 'G'
        elif c == 'G':
            rcread += 'C'
        elif c == 'T':
            rcread += 'A'
        else:
            rcread += 'A'
    return rcread[::-1]

def cdbg_read(f):
    idx, copy_num = 0, 0

    state = 0
    for line in ifilter(lambda x: len(x)>0, imap(string.strip, f)):
        if state == 0:
            m = re.match('>seq_(\d+)\s+(\d+)', line)
            if m:
                copy_num = int(m.group(2))
                state = 1
            else:
                raise Exception('invalid cdbg file')
        elif state == 1:
            yield idx, line, copy_num
            idx += 1
            state = 0

def component_read(f):
    idx, contigs, gaps = 0, [], []

    state = 0
    for line in imap(string.strip, f):
        if state == 0:
            m = re.match('>component\s+(\d+)', line)
            if m:
                state = 1
            else:
                raise Exception('invalid component file')
        elif state == 1:
            contigs = map(int, line.split())
            state = 2
        elif state == 2:
            gaps = map(int, line.split())
            yield idx, contigs, gaps
            idx += 1
            state = 0

USAGE = '%s [workspace] [cdbg] [omponent]'

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print USAGE % sys.argv[0]
        sys.exit(1)

    workspace, cdbg_file, component_file = sys.argv[1:4]

    cdbg_tbl, length_tbl, reverse_tbl = {}, {}, {} # i=>j

    with file(cdbg_file, 'r') as f:
        for idx, contig, copy_num in cdbg_read(f):
            #if copy_num > 1:
                #continue
            cdbg_tbl[contig] = idx
            length_tbl[idx] = len(contig)

    for contig, pos in cdbg_tbl.iteritems():
        if length_tbl[pos] < 100:
            continue
        if pos in reverse_tbl:
            continue
        rr = make_complement(contig)
        if rr in cdbg_tbl:
            reverse_tbl[pos] = cdbg_tbl[rr]
            reverse_tbl[cdbg_tbl[rr]] = pos

    print "[info] num of unique contig = %d" % (len(cdbg_tbl))
    print "[info] num of paired contig = %d" % (len(reverse_tbl))

    with file(os.path.join(workspace, 'misspaired.fa'), 'w') as f:
        for contig, pos in cdbg_tbl.iteritems():
            if pos not in reverse_tbl:
                print>>f, '>ctg'
                print>>f, contig 
                print>>f, make_complement(contig)

    scflens = {}
    sumL = 0
    scfs = []	# [[ids]] ids of ctgs
    ctg2scf = {}	#id of ctg 2 id of scf

    with file(component_file, 'r') as f:
        try:
            for idx, contigs, gaps in component_read(f):
                scfs.append((contigs, gaps))
                for contig in contigs:
                    ctg2scf[contig] = idx
                l = sum(length_tbl[contig] for contig in contigs) + sum(gaps)
                sumL += l
                scflens[idx] = l
        except:
            info = sys.exc_info()
            print info

    print "[info] num of scfs = %d" % (len(scfs))

    sortedLens = sorted(scflens.iteritems(), key=lambda x: x[1], reverse=True)
    with file(os.path.join(workspace, 'scflens_predict_0.data'), 'w') as f:
        for _, l in sortedLens:
            print>>f, l

    print "[info] sumG = %d" % (sumL)
    tmpl = 0
    for _, l in sortedLens:
        tmpl += l
        if tmpl >= sumL/2:
            print "[info] N50 = %d" % (l)
            break

    flag = [{} for i in scfs]
    for idx, _ in sortedLens:
        contigs, _ = scfs[idx]
        for contig in contigs:
            if not contig in reverse_tbl:
                continue
            rr = reverse_tbl[contig]
            if rr in ctg2scf:
                sid = ctg2scf[rr]
                if sid in flag[idx]:
                    flag[idx][sid] += 1
                else:
                    flag[idx][sid] = 1

    duplication = 0
    sumL = 0
    newscfs = []
    dflag = [0]*len(scfs)

    for idx, l in sortedLens:
        if dflag[idx] == 1:
            continue
        newscfs.append(idx)
        sumL += l
        for sid, w in flag[idx].iteritems():
            if w <= 1:
                continue
            if dflag[sid] == 1:
                duplication += 1
            dflag[sid] = 1

    print "[info] duplication = %d" % (duplication)
    print "[info] remain scfs = %d" % (len(newscfs))

    with file(os.path.join(workspace, 'component_last'), 'w') as f:
        count = 0
        for idx in newscfs:
            print>>f, '>component\t%d' % count
            contigs, gaps = scfs[idx]
            print>>f, ' '.join('%d' % i for i in contigs)
            print>>f, ' '.join('%d' % i for i in gaps)
            count += 1

    with file(os.path.join(workspace, 'scflens_predict_1.data'), 'w') as f:
        for idx in newscfs:
            print>>f, scflens[idx]

    print "[info] sumG = %d" % sumL
    tmpl = 0
    for index in newscfs:
        tmpl += scflens[index]
        if tmpl >= sumL/2:
            print "[info] N50 = %d" % scflens[index]
            break

