#!/bin/bash
#
# created by wangbing@ict.ac.cn
#
#################################################
if [ $# -lt 1 ]; then
    echo "Usage $0 [cfg.txt] ..."
    exit 1
fi

cwd=$(dirname `readlink -f $0`)
cd ${cwd}

ulimit -c unlimited

pyarcs="../src/ARCS.py"
for cfg_file in $@; do
    echo ${cfg_file}

    while read line
    do
        echo $line
        dbname=`echo $line | cut -d " " -f1`
        cfg=`echo $line | cut -d " " -f2`
        start_k=`echo $line | cut -d " " -f3`
        end_k=`echo $line | cut -d " " -f4`

        result_path="data/$dbname"
        mkdir -p $result_path

        for K in  `seq $start_k 2 $end_k` 
        do
            mkdir -p $result_path/$K"mer"
            $pyarcs -K $K -d $result_path/$K"mer"/  -s $cfg -e $K -O 200 -E 1 | tee $result_path/$K"mer/"$K"mer.log"
        done
    done < ${cfg_file}
done

