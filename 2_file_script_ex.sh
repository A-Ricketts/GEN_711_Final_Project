#! /bin/bash

source activate genomics

data=$1

for item in $data*R1*
 do
  reverse=$(ls $data* | grep -A1 "$item" | tail -n1)
  echo $item
  echo $reverse
done
