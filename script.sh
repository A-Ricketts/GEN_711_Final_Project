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

for item in $data*
 do
  fastqc $item -o $rawreads
done

sleep 5

#adpater and quality trimming

for item in $data*R1*
 do
  reverse=$(ls $data* | grep -A1 "$item" | tail -n1)
  trim_scriptV2.sh $item $reverse
done

sleep 5

#genome assembly

mkdir $data../spades_assembly

spadesout="$data../spades_assembly"
trimmed_reads="$data../trimmed-reads/"

for item in $(ls $trimmed_reads)
 do
if [[ "$item" == *R1* ]]
then
if [[ "$item" != *unpaired* ]]
then
  reverse=$(ls $trimmed_reads | grep -A1 "$item" | head -n2 | tail -n1)
  unpairedf=$(ls $trimmed_reads | grep "$item" | tail -n1)
  unpairedb=$(ls $trimmed_reads | grep -A1 "$unpairedf" | tail -n1)
  spades.py -1 $trimmed_reads$item -2 $trimmed_reads$reverse -s $trimmed_reads$unpairedf -s $trimmed_reads$unpairedb -o $spadesout -t 24
fi
fi
done

#remove unnecessary genome files

for item in $spadesout*
 do
  if [[ "$item" != *contigs.fasta ]]
   then
    if [[ "$item" != *spades.log ]]
     then
      rm $item
    fi
  fi
done

#genome assessment

mkdir $data../quast_results

$quast="$data../quast_results"

for item in $spadesout*
 do
  if [[ "$item" == *contigs.fasta ]]
   then
    quast.py $item -o $quast
  fi
done

#busco

mkdir $data../busco_results

$buscoo="$data../busco_results"

for item in $spadesout*
 do
  if [[ "$item" == *contigs.fasta ]]
   then
    busco -i $item -m genome -o $buscoo -l bacteria
  fi
done

#genome annotation

mkdir $data../prokka_output

$prokka="$data../prokka_output/"

for item in $spadesout*
 do
  if [[ "$item" == *contigs.fasta ]]
   then
    prokka $item --outdir $prokka --cpus 24 --mincontiglen 200
  fi
done

mkdir $data../protein_abundances

protein_ab="$data../protein_abundances"

for item in $prokka
 do
  if [[ "$item" == *PROKKA_*.gff ]]
   then
    grep -o "product=.*" $item | sed 's/product=//g' | sort | uniq -c | sort -nr > "$item"protein_abundances.txt
    mv "$item"protein_abundances.txt $protein_ab
  fi
done

###organism identification

#16S sequence

#BLAST


