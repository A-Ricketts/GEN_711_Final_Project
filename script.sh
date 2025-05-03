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
  sleep 5
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
      sleep 5
    fi
  fi
done

sleep 5

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

for item in $spadesout*contigs.fasta
 do
  quast.py $item -o $quast
  sleep 5
done

sleep 5

#busco

mkdir $data../busco_results

$buscoo="$data../busco_results"

for item in $spadesout*contigs.fasta
 do
  busco -i $item -m genome -o $buscoo -l bacteria
  sleep 5
done

sleep 5

#genome annotation

mkdir $data../prokka_output

$prokka="$data../prokka_output/"

for item in $spadesout*contigs.fasta
 do
  prokka $item --outdir $prokka --cpus 24 --mincontiglen 200
  sleep 5
done

sleep 5

mkdir $data../protein_abundances

protein_ab="$data../protein_abundances"

for item in $prokka*.gff
 do
  grep -o "product=.*" $item | sed 's/product=//g' | sort | uniq -c | sort -nr > "$item"protein_abundances.txt
  mv "$item"protein_abundances.txt $protein_ab
  sleep 5
done

sleep 5

###organism identification

#16S sequence

mkdir $data../16S_sequences

sixtS_seq="$data../16S_sequences"

for item in $prokka*.ffn
 do
  extract_sequences "16S ribosomal RNA" $item > "$item"16S_sequence.fasta
  mv "$item"16S_sequence.fasta $sixtS_seq
  sleep 5
done

sleep 5

#BLAST

mkdir $data../contigs_db

con_db="$data../contigs_db"

for item in $spadesout*contigs.fasta
 do
  makeblastdb -in $item -dbtype nucl -out $con_db
  sleep 5
done

mkdir $data../vs_contig

vs_contig="$data../vs_contig"

for item in $sixtS_seq
 do
  blastn -query $item -db $con_db -out "$item"16S_vs_contigs_6.tsv -outfmt 6
  mv "$item"16S_vs_contigs_6.tsv $vs_contigs
  sleep 5
done

for item in $spadesout*contigs.fasta
 do
  blob_blast.sh $item
done

#read mapping

###unsure of what this is

mkdir $data../raw_mapped

raw_map="$data../raw_mapped"
reads=$(ls $trimmed_reads)
a=1

for item in $spadesout*contigs.fasta
 do
  forward=$(grep "R1" $reads | head -n"$a" | tail -n1)
  reverse=$(grep "R2" $reads | head -n"$a" | tail -n1)
  bwa index $item
  bwa mem -t 24 $item $forward $reverse > "$item"raw_mapped.sam
  mv "$item"raw_mapped.sam $raw_map
  a=$($a - 1)
  sleep 5
done

mkdir $data../sorted_mapped

sort_map="$data../sorted_mapped"

for item in $raw_map
 do
  samtools view -@ 24 -Sb $item | samtools sort -@ 24 -o "$item"sorted_mapped.bam
  mv "$item"sorted_mapped.bam $sort_map
  sleep 5
done

mkdir $data../bed_cov_out

bed_out="$data../bed_cov_out"

for item in $sort_map
 do
  samtools flagstat $item
  samtools index $item
  bedtools genomecov -ibam $item > "$item"coverage.out
  mv "$item"coverage.out $bed_out
  sleep 5
done

































