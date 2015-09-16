#!/bin/bash
#write by wangbing
#devide LP to small pieces and run glpk on the small LPs
#
PWD=$(dirname $(readlink -f $0))
if [ $# -ne 3 ]; then
	echo "usage: $0 bigLPfile outSmallLPPath resultOutPath" 
	exit 1
fi

if [ ! -d $2 ]; then
	mkdir -p $2
fi

if [ ! -d $3 ]; then
	mkdir -p $3
fi

${PWD}/arcs divideLP -i $1 -d $2 

for f in $(ls $2); do
	if [ `echo $f | grep "math" | wc -l` -eq 0 ]; then
		continue
	fi
	name=`echo $f | sed "s/math/out/g"`
#	echo "glpsol --math $2/$f -o $3/$name"
	glpsol --math $2/$f -o $3/$name > /dev/null
done
