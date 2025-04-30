#! /bin/bash

source activate genomics

#estimate of coverage of reads

data="~/Final_Project/Data/"

for item in ~/Final_Project/Data/*R1*
  do
   num_reads=$(zgrep -c "^@" $item)
   bp=$(($num_reads * 250 * 2 / 7000000))
   if [ $bp -ge 70 ]
    then
     echo $item has good coverage
    else
     echo $item does not have good coverage
   fi
done

sleep 5

#read quality

mkdir ~/Final_Project/fastqc_raw-reads

rawreads="~/Final_Project/fastqc_raw-reads"

for item in $data
 do
  fastqc $item -o $rawreads
done

sleep 5

#adpater and quality trimming

mkdir ~/Final_Project/trimmed_reads

trimreads="~/Final_Project/trimmed_reads"
trimqc="~/Final_Project/trimmed_reads_fastqc"

for item in $data
 do
  trim_file=$(trim_scriptV2.sh $item)
  mv $trim_file $trimreads
  fastqc $trim_file -o $trimqc
done

sleep 5

#genome assembly

