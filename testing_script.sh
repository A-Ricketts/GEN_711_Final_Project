#! /bin/bash
source activate genomics

spadesout=$1

fasta=$(ls $spadesout | grep  | head -n1 | tail -n1)
echo $fasta
