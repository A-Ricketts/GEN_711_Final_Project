#! /bin/bash

source activate genomics

data=$1
hey="hey"
spadesout="$data../spades_assembly/"

for item in $spadesout*
 do
if [[ "$item" == *contigs.fasta ]]
then
  echo $item
fi
done

for item in $spadesout*contigs.fasta
do 
echo $item
done
