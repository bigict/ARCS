#!/bin/sh
#!/bin/bash
echo "enter shell!!!"
echo "begin map reads to training data..."

echo `date` > time
:<<BLOCK
i=1
while getopts "1:2:t:g:" arg
do 
	case $arg in 
		1)
		read_file_name1=$OPTARG
		;;
		2)
		read_file_name2=$OPTARG
		;;
		t)
		training=$OPTARG
		;;
		g)
		multi_gap_seq_to_PerM=$OPTARG
		;;
		?)
		echo "Invalid parameter"
		exit 1
		;;
	esac
done


mkdir ARCS_TEMP
cp $training ./ARCS_TEMP/
cp $multi_gap_seq_to_PerM ./ARCS_TEMP/

dir="./ARCS_TEMP/"
file=$multi_gap_seq_to_PerM
training=$training

echo "map uniqEdgeForTraining" >> time
echo `date` >> time

./perm10K $dir/$training -v 10 -a  -1 $read_file_name1 -2 $read_file_name2 -s $dir/$training.index --outputFormat sam -o $dir/$training.sam


echo "map UniqEdgeForTraining end" >> time
echo `date` >> time

echo "map multi_gap_seq_to_PerM" >> time
echo `date` >> time

./perm10K $dir/$file -v 10 -a  -1 $read_file_name1 -2 $read_file_name2 -s $dir/$file --outputFormat sam -o $dir/$file.sam

echo "map multi_gap_seq_to_PerM end" >> time
echo `date` >> time
BLOCK
