#!/bin/bash


#./ARCS.py -s ../../../../E.coli.cfg -d ../../../result/Ecoli/29mer/ -K 29 -p 5 > ../../../result/Ecoli/29mer/29mer.log &
#
#./ARCS.py -s ../../../../E.coli.cfg -d ../../../result/Ecoli/27mer/ -K 27 -p 2 > ../../../result/Ecoli/27mer/27mer.log &
#
#./ARCS.py -s ../../../../E.coli.cfg -d ../../../result/Ecoli/25mer/ -K 25 -p 2 > ../../../result/Ecoli/25mer/25mer.log &
#
#./ARCS.py -s ../../../../E.coli.cfg -d ../../../result/Ecoli/ -K 31 -p 2 > ../../../result/Ecoli/31mer.log &


for k in `seq 21 2 35`;do

	mkdir ../../../result/sta/${k}mer
	./ARCS.py -s ../../../../Sta.cfg -d ../../../result/sta/${k}'mer'/ -K $k -p 2 > ../../../result/sta/${k}'mer'/${k}'mer'.log  &
done
