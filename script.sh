#! /bin/bash

conda activate genomics

#estimate of coverage of reads

data="~/Final_Project/Data/"

##something is wrong!!
for item in ~/Final_Project/Data/*R1*
  do
   num_reads=$(zgrep -c "^@" $item)
   bp=$($num_reads * 250 * 2 / 7000000)
   if [ $bp -ge 70 ]
    then
     echo $item has good coverage
    else
     echo $item does not have good coverage
   fi
done

#read quality

mkdir ~/Final_Project/fastqc_raw-reads

rawreads="~/Final_Project/fastqc_raw-reads"

for item in $data
 do
  fastqc $item $rawreads
done

sleep 5

#adpater and quality trimming

mkdir ~/Final_Project/trimmed_reads

trimreads="~/Final_Project/trimmed_reads"


