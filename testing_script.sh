#! /bin/bash

source activate genomics

#estimate of coverage of reads

data="~/Final_Project/Data/"

##something is wrong!!
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

