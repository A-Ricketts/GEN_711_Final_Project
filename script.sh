#! /bin/bash

source activate genomics

#estimate of coverage of reads

data=$1

touch $data../coverage.txt
coverage="$data../coverage.txt"

for item in $data*R1*
  do
   num_reads=$(zgrep -c "^@" $item)
   bp=$(($num_reads * 250 * 2 / 7000000))
   if [ $bp -ge 70 ]
    then
     echo $item has good coverage >> $coverage
    else
     echo $item does not have good coverage >> $coverage
   fi
done

sleep 5

#read quality

mkdir $data../fastqc_raw-reads

rawreads="$data../fastqc_raw-reads"

for item in $data
 do
  fastqc $item -o $rawreads
done

sleep 5

#adpater and quality trimming

mkdir $data../trimmed_reads

trimreads="$data../trimmed_reads"

for item in $data*R1*
 do
  reverse=$(ls $data* | grep -A1 "$item" | tail -n1)
  trim_scriptV2.sh $item $reverse
  mv *paired* $trimreads
done

sleep 5

#genome assembly

