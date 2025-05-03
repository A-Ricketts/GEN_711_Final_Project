#! /bin/bash

source activate genomics

data=$1

for item in $data*R1*
  do
   echo $item
done
