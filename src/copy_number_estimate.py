#!/usr/bin/python
# -*- coding: utf-8 -*-

import os,sys
from itertools import ifilter,imap
import string

def Worker_read(f):
    state = 0
    for line in ifilter(lambda x: len(x)>0, imap(string.strip, sys.stdin)):
        try:
            if state == 0:
                i, j, seq = line.split('\t')[:3]
                state = 1
            elif state == 1:
                coverage = map(float, line.split())
                if len(coverage) > 0:
                    coverage = sum(coverage) / len(coverage)
                else:
                    coverage = 0
                yield int(i), int(j), seq, coverage
                state = 0
        except:
            print>>sys.stderr, 'invalid line: %s' % (line)

def Worker_run(config, contig, component):
    config_file, component_file = file(os.path.join(config['workdir'], contig), 'wb'), file(os.path.join(config['workdir'], component), 'wb')

    contig_index, component_index, edge_index = 0, 0, 0
    for _, _, seq, coverage in Worker_read(sys.stdin):
        copy_num = int(round(coverage / float(config['lambda'])))
	if copy_num >= 0:
            print>>config_file, '>seq_%d\t%d' % (contig_index, copy_num)
            print>>config_file, seq

            if copy_num <= 1:
	        print>>component_file, '>component\t%d' % (component_index)
                print>>component_file, '%d' % (contig_index)
                print>>component_file, ''
                component_index += 1

            contig_index += 1
        edge_index += 1

    print '\tNumber of condensed edges = %d' % (edge_index)
    print '\tNumber of Contigs = %d' % (contig_index)
    print '\tNumber of Uniques = %d' % (component_index)

USAGE = '%s -p [contig_parameter_file] [output_fasta_file] [component0]'

if __name__ == '__main__':
    import getopt

    config = {'workdir': '.'}
    opts, args = getopt.getopt(sys.argv[1:], "i:p:d:h")
    for o, a in opts:
        if o == '-h':
            print USAGE % (sys.argv[0])
            sys.exit(1)
        elif o == '-d':
            config['workdir'] = a
        elif o == '-p':
            with file(a) as f:
                for line in ifilter(lambda x: len(x) > 0 and not x.startswith('#'), imap(string.strip, f)):
                    k, v = imap(string.strip, line.split('=', 1))
                    config[k] = v

    if len(args) == 2:
        Worker_run(config, *args)
    else:
        print USAGE % (sys.argv[0])
