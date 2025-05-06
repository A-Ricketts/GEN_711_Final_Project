#! /bin/bash
source activate genomics

data=$1
spadesout=$1

fasta=$(head -n1 $spadesout*contigs.fasta | tail -n1)
echo $fasta
