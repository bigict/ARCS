#! /usr/bin/python

Usage = '''
Author : Zheng QuanGang
Function : 

$python dbfindDFS.py	DFS/BFS/DFSEXT
                        <condensed de Bruijn graph> 
			<cdbg>
			<workspace>
                        <start contig> <end contig>
                        //
'''

import sys

mode = "DFS"
if len(sys.argv) != 7:
	print Usage
	sys.exit()
else:
	debruijn = sys.argv[2]
	cdbg_file = sys.argv[3]
	prefix = sys.argv[4]
	s1 = int(sys.argv[5])
	s2 = int(sys.argv[6])
        mode = sys.argv[1]
	
def reversecomplement( read ):
	rcread = ''
	for i in range(len(read)):
		if read[len(read)- i - 1] == 'A':
			rcread = rcread + 'T'
			continue
		if read[len(read)- i - 1] == 'T':
			rcread = rcread + 'A'
			continue
		if read[len(read)- i - 1] == 'C':
			rcread = rcread + 'G'
			continue
		if read[len(read) - i - 1] == 'G':
			rcread = rcread + 'C'
			continue
		if read[len(read) - i - 1] == 'N':
			rcread = rcread + 'N'
        return rcread

de = {}
rcontigs = []
contigs = []
ref_pos = {}
ctg_pos = []
node_index = {}
node_pre_index = {}
cov = {}
cov_avg = {}


def readcdbg():
	global rcontigs
	global cdbg_file
	rcontigs = []
	cin = open(cdbg_file,"r")
	try:
		while 1:
			line =cin.next().strip()
			line = cin.next().strip()
			rcontigs.append(line)
	except:
		cin.close()

def readdeb():
	global contigs
	global rcontigs
	global debruijn
	global de
        global cov_avg
	contigs = []
	de = {}
	index = 0
	cindex = 0 
	cin = open(debruijn,"r")
	try:
	    while 1:
	        line = cin.next().strip()
	        array = line.split()
	        i = int(array[0])
	        j = int(array[1])
	        contigs.append(array[2])
	        line1 = cin.next().strip()
	        array1 = line1.split()
	        if array[2] in rcontigs:
	            de[(i,j)] = index
                    cov[index] = array1
	            index += 1
	            cindex += 1
	        else:
                    de[(i,j)] = "x"+str(cindex)
                    cov["x"+str(cindex)] = array1
	            cindex += 1
	except:
	    info = sys.exc_info()
	    print info[0],info[1]
	    cin.close()
	print "Step 1\tNumber of contigs = "+str(index)
	print "      \tNumber of edges = "+str(cindex)
        for i in cov:
            sum1 = 0.0
            for j in cov[i]:
                sum1 += int(j)
            cov_avg[i] = int(sum1 / len(cov[i]))



def debindex():
	global node_index
        global node_pre_index
	global de
	index = 0
	node_index = {}
	for i,j in de:
	    if i in node_index:
	        node_index[i].append(j)
	    else:
	        node_index[i] = []
	        node_index[i].append(j)
	    if j in node_pre_index:
	        node_pre_index[j].append(i)
	    else:
	        node_pre_index[j] = []
	        node_pre_index[j].append(i)

de1node = {}
de2node = {}
def edgeindex():
	global de1node
	global de2node
	global de
	de1node = {}
	de2node = {}
	for i,j in de:
	    if type(de[(i,j)]) == int:
	        de1node[int(de[(i,j)])] = i
	        de2node[int(de[(i,j)])] = j

pathseq = []
paths = []
path1 = []
path2 = []
cov0 = []

def DFS(a,b,c,ptmp,path,step):
    global node_index
    #print step 
    if c == b:                     
        print "found one path"
        path.append(ptmp[:])          
        return
    if step >= 100:
        return
    if c in node_index:
        ptmp.append(c)
        for k in node_index[c]:             
            DFS(a,b,k,ptmp,path,step+1)
        tmp = ptmp.pop()
    else:
        return

def DFSfind():
	global s1
	global s2
	global paths 
        global cov0
        global cov
        paths = []
	start = de2node[s1]
	end = de1node[s2]
        print de1node[s1],de2node[s1],s1,len(rcontigs[s1])
        print de1node[s2],de2node[s2],s2,len(rcontigs[s2])
	cin = open(prefix+"/Cyto_s",'w')
	fout = open(prefix+"/"+str(s1)+"_"+str(s2)+"gap.fa","w")
        cin.write("i\tj\tcontig_id|len\n")
	cin.write(str(de1node[s1])+"\t"+str(de2node[s1])+"\t"+str(s1)+"|"+str(len(rcontigs[s1]))+"\n")
	cin.write(str(de1node[s2])+"\t"+str(de2node[s2])+"\t"+str(s2)+"|"+str(len(rcontigs[s2]))+"\n")
	ptmp = []
	k = start
        step = 0
	DFS(start,end,k,ptmp,paths,step)
        print "Paths = ",len(paths)
        ci = 0
        seq = ""
	for i in paths:
                covs = []
                print "path"+str(ci),i,len(i)
                fout.write(">seq_0_"+str(ci)+"\n")
                ci += 1
                fout.write(rcontigs[s1])
                seq = rcontigs[s1]
                covs += cov[de[(de1node[s1],de2node[s1])]]
		for j in range(len(i)-1):
			x = i[j]
			y = i[j+1]
                        covs += cov[de[(x,y)]]
	                if type(de[(x,y)]) == int:
	                    cin.write(str(x)+"\t"+str(y)+"\t"+str(de[(x,y)])+"|"+str(len(rcontigs[de[(x,y)]]))+"|"+str(cov_avg[de[(x,y)]])+"\n")
                            fout.write(rcontigs[de[(x,y)]][30:])
                            seq += rcontigs[de[(x,y)]][30:]
	                else:
	                    index = int(de[(x,y)][1:])
	                    cin.write(str(x)+"\t"+str(y)+"\t"+str(de[(x,y)])+"|"+str(len(contigs[index]))+"|"+str(cov_avg[de[(x,y)]])+"\n")
                            fout.write(contigs[index][30:])
                            seq += contigs[index][30:]
                if len(i) == 0:
                    continue
		x = i[len(i)-1]
		y = end
                covs += cov[de[(x,y)]]
	        if type(de[(x,y)]) == int:
	            cin.write(str(x)+"\t"+str(y)+"\t"+str(de[(x,y)])+"|"+str(len(rcontigs[de[(x,y)]]))+"\n")
                    fout.write(rcontigs[de[(x,y)]][30:])
                    seq += rcontigs[de[(x,y)]][30:]
	        else:
	            index = int(de[(x,y)][1:])
	            cin.write(str(x)+"\t"+str(y)+"\t"+str(de[(x,y)])+"|"+str(len(contigs[index]))+"\n")
                    fout.write(contigs[de[(x,y)]][30:])
                    seq += contigs[index][30:]
                fout.write(rcontigs[s2][30:])
                covs += cov[de[(de1node[s2],de2node[s2])]]
                seq += rcontigs[s2][30:]
                cov0.append(covs)
                fout.write("\n")
                pathseq.append(seq)
                #print seq
	cin.close()
        fout.close()

def BFSfind():
	global prefix
	global de1node
	global de2node
	global rcontigs
	global contigs
	global node_index
	global de
	global s1
	global s2
	cin = open(prefix+"/Cyto_s",'w')
	stack = []
	start = de2node[s1]
	end = de1node[s2]
	print de1node[s1],de2node[s1],s1,len(rcontigs[s1])
	print de1node[s2],de2node[s2],s2,len(rcontigs[s2])
	cin.write("i\tj\tcontig_id|len\n")
	cin.write(str(de1node[s1])+"\t"+str(de2node[s1])+"\t"+str(s1)+"|"+str(len(rcontigs[s1]))+"\n")
	cin.write(str(de1node[s2])+"\t"+str(de2node[s2])+"\t"+str(s2)+"|"+str(len(rcontigs[s2]))+"\n")
	stack.append(start)
	depth = 0
	count = 1
	tmp = 0
	while len(stack) != 0:
	    k = stack.pop(0)
	    count -= 1
	    if k == end:
	        print "found"
	        if count <= 0:
	            count = tmp
	            tmp = 0
	            depth += 1
	        #continue
	    if depth >= 30:
	        continue
	    if k in node_index:
	        for z in node_index[k]:
	            if type(de[(k,z)]) == int:
	                cin.write(str(k)+"\t"+str(z)+"\t"+str(de[(k,z)])+"|"+str(len(rcontigs[de[(k,z)]]))+"|"+str(cov_avg[de[(k,z)]])+"\n")
	            else:
	                index = int(de[(k,z)][1:])
	                cin.write(str(k)+"\t"+str(z)+"\t"+str(de[(k,z)])+"|"+str(len(contigs[index]))+"|"+str(cov_avg[de[(k,z)]])+"\n")
	
	            stack.append(z)
	            tmp += 1
	    if count <= 0:
	        count = tmp
	        tmp = 0
	        depth += 1
	cin.close()
    
def DFS_backward(c,ptmp,path,lens):
    #global node_index
    global node_pre_index
    index = node_pre_index
    if c in index:
        for k in index[c]:
            cont = ""
            ptmp.append(k)
            if type(de[(k,c)]) == int:
                cont = rcontigs[de[(k,c)]]
            else:
                cont = contigs[int(de[(k,c)][1:])]
            lens = lens + len(cont) - 30
            if lens >= 100:
                path.append(ptmp[0:])
                #print ptmp[:]
                #continue
            else:
                DFS_backward(k,ptmp,path,lens)
            tmp = ptmp.pop()
            lens = lens - len(cont) + 30
            
def DFS_forward(c,ptmp,path,lens):
    global node_index
    index = node_index
    if c in index:
        for k in index[c]:
            cont = ""
            ptmp.append(k)
            if type(de[(c,k)]) == int:
                cont = rcontigs[de[(c,k)]]
            else:
                cont = contigs[int(de[(c,k)][1:])]
            lens = lens + len(cont) - 30
            if lens >= 100:
                path.append(ptmp[0:])
                #continue
            else:
                DFS_forward(k,ptmp,path,lens)
            tmp = ptmp.pop()
            lens = lens - len(cont) + 30


longseq = []
def DFSfind_ext():
	global s1
	global s2
        global node_index
        global node_pre_index
        global path1
        global path2
        global pathseq
        global paths
        global cov
        global longseq
        longseq = []
	a = de1node[s1]
	c = de2node[s2]
        #print de1node[s1],de2node[s1],s1,len(rcontigs[s1])
        #print de1node[s2],de2node[s2],s2,len(rcontigs[s2])
	#cin = open(prefix+"/Cyto_s",'w')
        #cin.write("i\tj\tcontig_id|len\n")
	#cin.write(str(de1node[s1])+"\t"+str(de2node[s1])+"\t"+str(s1)+"|"+str(len(rcontigs[s1]))+"\n")
	#cin.write(str(de1node[s2])+"\t"+str(de2node[s2])+"\t"+str(s2)+"|"+str(len(rcontigs[s2]))+"\n")
	path1 = []
        seq1 = []
        cov1 = []
        path2 = []
        seq2 = []
        cov2 = []
	ptmp = [a]
	
        DFS_backward(a,ptmp,path1,0)
        print "Path1 = ",len(path1)
        
        seq = ""
        covs = []
        ci = 0
	for p in path1:
                seq = ""
                covs = []
                lis = [len(p)-i-1 for i in range(len(p))]
                if len(lis) <= 1:
                    continue
                x = p[lis[0]]
                y = p[lis[1]]
                #print x,y
                covs += cov[de[(x,y)]]
	        if type(de[(x,y)]) == int:
                    seq += rcontigs[de[(x,y)]]
	        else:
	            index = int(de[(x,y)][1:])
                    seq += contigs[index]
                for j in range(len(lis)-2):
			x = p[lis[j+1]]
			y = p[lis[j+2]]
                        #print x,y
                        covs += cov[de[(x,y)]]
	                if type(de[(x,y)]) == int:
                            seq += rcontigs[de[(x,y)]][30:]
	                else:
	                    index = int(de[(x,y)][1:])
                            seq += contigs[index][30:]

                #print y
                seq1.append(seq)
                #print seq
                cov1.append(covs)
                ci += 1

        ptmp = [c]
        DFS_forward(c,ptmp,path2,0)
        print "Path2 = ",len(path2)
        seq = ""
        ci = 0
	for p in path2:
                seq = ""
                covs = []
                if len(p) <= 1:
                    continue
                x = p[0]
                y = p[1]
                covs += cov[de[(x,y)]]
	        if type(de[(x,y)]) == int:
                    seq += rcontigs[de[(x,y)]]
	        else:
	            index = int(de[(x,y)][1:])
                    seq += contigs[index]
                for j in range(len(p)-2):
			x = p[j+1]
			y = p[j+2]
                        covs += cov[de[(x,y)]]
	                if type(de[(x,y)]) == int:
                            seq += rcontigs[de[(x,y)]][30:]
	                else:
	                    index = int(de[(x,y)][1:])
                            seq += contigs[index][30:]
                seq2.append(seq)    
                cov2.append(covs)
                ci += 1
        
	fout = open(prefix+"/"+str(s1)+"_"+str(s2)+"longseq.fa","w")
        fout1 = open(prefix+"/"+str(s1)+"_"+str(s2)+"longcov.fa","w")
        lseq = ''
        lcov = []
        lseqs = []
        lcovs = []
        ci = 0
        #for i in seq1:
        #    print i
        for p in range(len(pathseq)):
            lcov = cov0[p]
            lseq = pathseq[p]
            print len(lcov)
            if len(seq1) > 0:
                for p1 in range(len(seq1)):
                    lseq = seq1[p1][-130:] + pathseq[p][30:]
                    #print seq1[p1][-30:],pathseq[p][:30]
                    lcov = cov1[p1][-100:] + cov0[p]
                    if len(seq2) > 0:
                        for p2 in range(len(seq2)):
                            ls = lseq + seq2[p2][30:130]
                            lc = lcov + cov2[p2][:100]
                            lseqs.append(ls)
                            lcovs.append(lc)
                            fout.write(">seq"+str(p1)+str(p)+str(p2)+"_0_"+str(ci)+"\n"+ls+"\n")
                            fout1.write(">seq"+str(p1)+str(p)+str(p2)+"_0_"+str(ci)+"\n")
                            for i in lc:
                                fout1.write(str(i)+" ")
                            ci += 1
                            fout1.write("\n")
                    else:
                        lseqs.append(lseq)
                        lcovs.append(lcov)
                        fout.write(">seq"+str(p1)+str(p)+str('n')+"_0_"+str(ci)+"\n"+lseqs+"\n")
                        fout1.write(">seq"+str(p1)+str(p)+str('n')+"_0_"+str(ci)+"\n")
                        for i in lcovs:
                            fout1.write(str(i)+" ")
                        ci += 1
                        fout1.write("\n")
            else:
                if len(seq2) > 0:
                    for p2 in range(len(seq2)):
                        ls = lseq + seq2[p2][30:130]
                        lc = lcov + cov2[p2][:100]
                        lseqs.append(ls)
                        lcovs.append(lc)
                        fout.write(">seq"+str('n')+str(p)+str(p2)+"_0_"+str(ci)+"\n"+ls+"\n")
                        fout1.write(">seq"+str('n')+str(p)+str(p2)+"_0_"+str(ci)+"\n")
                        for i in lc:
                            fout1.write(str(i)+" ")
                        ci += 1
                        fout1.write("\n")
        fout.close()
        fout1.close()


readcdbg()
readdeb()
debindex()
edgeindex()
if mode == "DFS":
    DFSfind()
elif mode == "BFS":
    BFSfind()
elif mode == "DFSEXT":
    print "[MODE] dfs extend"
    DFSfind()
    DFSfind_ext()


