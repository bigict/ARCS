#!/usr/bin/env python3

import sys

def get_len(f):
    i = 0
    com_len = []
    com_cov = []
    for line in filter(lambda x:len(x)>0, map(str.strip, f)):
        if i % 2 == 0:
            a = line.split()
            com_cov.append(int(a[1]))
        if i % 2 == 1:
            com_len.append(len(line))
        i += 1
    return com_len, com_cov
            
def get_start_num(f):
    line_num = 0
    for line in f:
        line_num += 1
    assert line_num % 3 == 0
    return line_num / 3 + 1
if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit("Usage: %s component_0 cdbg.file component_2 cov_threshold" % sys.argv[0])

    with open(sys.argv[2], "r") as f:
        com_len, com_cov = get_len(f)
    with open(sys.argv[1], "r") as f:
        i = 0
        uniq = set()
        for line in f:
            if i % 3 == 1:
                uniq.add( int(line) )
            i += 1
    with open(sys.argv[3], "r") as f:
        start = get_start_num(f)

    with open(sys.argv[3], "a") as f:
        for i in range( len(com_len) ):
            if i not in uniq and ((com_len[i] > 100 and com_cov[i] <= int(sys.argv[4])) or com_len[i] > 200):
                f.write(">component %d\n" % start)
                f.write("%d\n" % i)
                f.write("\n")
                start += 1

