#!/usr/bin/python

import sys
import os
import re
from itertools import ifilter,imap
import string

def Worker_run(config):
    workspace, iteration, edge_num = config['workspace'], int(config['iteration']), int(config['EDGE_CLUSTER_NUM'])
    print 'edge num %d' % (edge_num)

    pos_array = [0] * edge_num

    root = os.path.join(workspace, 'smallLPResults_' + str(config['iteration']))
    for fname in os.listdir(root):
        row_num, column_num = 0, 0
        with file(os.path.join(root, fname), 'r') as f:
            for line in ifilter(lambda x: len(x)>0, imap(string.strip, f)):
                m = re.match('([a-zA-Z0-9_-]+):\s*(.*)', line)
                if m:
                    if m.group(1)   == 'Rows':    row_num    = int(m.group(2))
                    elif m.group(1) == 'Columns': column_num = int(m.group(2))
                else:
                    m = re.match('(\d)+\s+x_(\d+)\s+\w+\s+(\-*\d+).*', line)
                    if m:
                        pos_array[int(m.group(2))] = int(m.group(3))
    with file(os.path.join(workspace, 'edge_cluster_pos_%d' % iteration), 'w') as f:
        for pos in pos_array:
            print>>f, '%d' % (pos)

USAGE = '%s -p [scaffold_parameter_file] -d [workspace] -i [iteration]'

if __name__ == '__main__':
    import getopt

    config = {}
    opts, args = getopt.getopt(sys.argv[1:], "d:i:p:h")
    for o, a in opts:
        if o == '-h':
            print USAGE % (sys.argv[0])
            sys.exit(1)
        elif o == '-d':
            config['workspace'] = a
        elif o == '-i':
            config['iteration'] = int(a)
        elif o == '-p':
            with file(a) as f:
                for line in ifilter(lambda x: len(x) > 0 and not x.startswith('#'), imap(string.strip, f)):
                    k, v = imap(string.strip, line.split('=', 1))
                    config[k.strip()] = v.strip()

    if config.has_key('EDGE_CLUSTER_NUM') and config.has_key('workspace') and config.has_key('iteration'):
        Worker_run(config, *args)
    else:
        print USAGE % (sys.argv[0])
