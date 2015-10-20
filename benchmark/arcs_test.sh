#!/bin/bash
if [ $# != 1 ]; then
    echo "Usage $0 cfg.txt"
    exit 1
fi

all_cfg_file=$1
arcs="/home/wangbing/ARCS_v1/src/ARCS.py"
arcs="/home/wangbing/ARCS_V1.1/ARCS_new/src/ARCS.py"
arcs="/home/wangbing/ARCS_V1.2/src/ARCS_new.py"
while read line
do
    echo $line
    dbname=`echo $line | cut -d " " -f1`
    cfg=`echo $line | cut -d " " -f2`
    start_k=`echo $line | cut -d " " -f3`
    end_k=`echo $line | cut -d " " -f4`

    result_path="/home/wangbing/result/$dbname"
    mkdir -p $result_path

    ulimit -c unlimited
    for K in  `seq $start_k 2 $end_k` 
    do
        mkdir -p $result_path/$K"mer"
        $arcs -K $K -d $result_path/$K"mer"/  -s $cfg -e $K -O 200 -E 1 | tee $result_path/$K"mer/"$K"mer.log"
    done
done < $all_cfg_file
