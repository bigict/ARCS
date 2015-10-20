#!/bin/bash
if [ $# != 1 ]; then
    echo "Usage $0 cfg.txt"
    exit 1
fi

all_cfg_file=$1
arcs="/home/wangbing/ARCS_v1/src/ARCS.py"

while read line
do
#    echo $line
    dbname=`echo $line | cut -d " " -f1`
    cfg=`echo $line | cut -d " " -f2`
    start_k=`echo $line | cut -d " " -f3`
    end_k=`echo $line | cut -d " " -f4`

    result_path="/home/wangbing/result/$dbname"

    echo $dbname
    for K in  `seq $start_k 2 $end_k` 
    do
        dir=$result_path/$K"mer"
        max=`head -4 $dir/$K"mer".final_result | tail -1`
        n50=`head -24 $dir/$K"mer".final_result | tail -1 | cut -d ":" -f2`
        n90=`head -25 $dir/$K"mer".final_result | tail -1 | cut -d ":" -f2`
        time=`tail -1 $dir/$K"mer".log | cut -d " " -f 5`
        if [ -f /home/zheng/new_ARCS_gage_wangbing/$dbname/$K"mer"/out.report ]; then
            Gco=`head -12 /home/zheng/new_ARCS_gage_wangbing/$dbname/$K"mer"/out.report | tail -1 | cut -d "(" -f2 | cut -d "%" -f1`
            Sco=`head -12 /home/zheng/new_ARCS_gage_wangbing/$dbname/$K"mer"/out.report | tail -1 | cut -d "(" -f3 | cut -d "%" -f1`
            echo "\\hline"
            echo $K"mer" "&" $max "&" $n50 "&" $n90 "&" $Gco"\\%" "&" $Sco"\\%" "&" $time "\\\\"
        else
            echo "\\hline"
            echo $K"mer" "&" $max "&" $n50 "&" $n90 "&" $time "\\\\"
        fi
    done
done < $all_cfg_file
