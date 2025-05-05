#! /bin/bash
source activate genomics

data=$1

for item in $data*
do
echo $item
done
