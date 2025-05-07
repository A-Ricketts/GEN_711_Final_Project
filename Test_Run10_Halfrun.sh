#! /bin/bash

source activate genomics

data=$1
spadesout="$data../spades_assembly/"
sixtS_seq="$data../16S_sequences/"
trimmed_reads="$data../trimmed-reads/"

#BLAST

con_db="$data../contigs_db/"

mkdir $data../vs_contig

vs_contig="$data../vs_contig/"
d=1

for item in $sixtS_seq*
 do
  cont_dba=$(ls $con_db* | grep "contigs_db.nsq" | head -n"$d" | tail -n1)
  blastn -query $item -db $cont_dba -out 16S_vs_contigs_6.tsv -outfmt 6
  mv 16S_vs_contigs_6.tsv "$item"_16S_vs_contigs_6.tsv
  mv "$item"_16S_vs_contigs_6.tsv $vs_contigs
  d=$(($d + 1))
  sleep 5
done
sleep 5

mkdir $data../megablast_out

mega_out="$data../megablast_out/"

for item in $spadesout*contigs.fasta
 do
  blob_blast.sh $item
  mv $data*megablast.out $mega_out
done

sleep 5
#read mapping

mkdir $data../raw_mapped
raw_map="$data../raw_mapped/"
#reads=$(ls $trimmed_reads*)
a=1

for item in $spadesout*contigs.fasta
 do
  forward=$(grep "R1" $trimmed_reads* | grep -v "unpaired" | head -n"$a" | tail -n1)
  reverse=$(grep "R2" $trimmed_reads* | grep -v "unpaired" | head -n"$a" | tail -n1)
  bwa index $item
  bwa mem -t 24 $item $forward $reverse > raw_mapped.sam
  j=1
  fasta=$(ls $spadesout*contigs.fasta | head -n"$j" | tail -n1)
  mv raw_mapped.sam "$fasta"_raw_mapped.sam
  mv "$fasta"_raw_mapped.sam $raw_map
  a=$(($a + 1))
  j=$(($j + 1))
  sleep 5
done

sleep 5

mkdir $data../sorted_mapped

sort_map="$data../sorted_mapped/"

for item in $raw_map*
 do
  samtools view -@ 24 -Sb $item | samtools sort -@ 24 -o sorted_mapped.bam
  i=1
  fasta=$(ls $spadesout*contigs.fasta | head -n"$i" | tail -n1)
  mv sorted_mapped.bam "$fasta"_sorted_mapped.bam
  mv "$fasta"_sorted_mapped.bam $sort_map
  i=$(($i + 1))
  sleep 5
done

sleep 5
mkdir $data../flagstat
mkdir $data../bed_cov_out

flag="$data../flagstat/"
bed_out="$data../bed_cov_out/"

for item in $sort_map*
 do
  samtools flagstat $item > "$item"_flagstat.txt
  mv "$item"_flagstat.txt $flag
  samtools index $item
  bedtools genomecov -ibam $item > coverage.out
  h=1
  fasta=$(ls $spadesout*contigs.fasta | head -n"$h" | tail -n1)
  mv coverage.out "$fasta"_coverage.out
  mv "$fasta"_coverage.out $bed_out
  h=$(($h + 1))
  sleep 5
done

sleep 5

mkdir $data../cov_table

cov_table="$data../cov_table/"
b=1

for item in $bed_out*
 do
  fasta=$(ls $spadesout*contigs.fasta | head -n"$b" | tail -n1)
  gen_input_table.py --isbedfiles $fasta $item > coverage_table.tsv
  mv coverage_table.tsv "$fasta"_coverage_table.tsv
  mv "$fasta"_coverage_table.tsv $cov_table
  b=$(($b + 1))
  sleep 5
done

sleep 5

#non-target contig removal

mkdir $data../blob_out

blob_out="$data../blob_out/"
c=1

for item in $spadesout*contigs.fasta
 do
  sort_bam=$(ls $sort_map* | head -n"$c" | tail -n1)
  mega=$(ls $mega_out* | head -n"$c" | tail -n1)
  blobtools create -i $item -b $sort_bam -t $mega -o blob_out
  mv blob_out* $blob_out
  c=$(($c + 1))
  sleep 5
done

sleep 5

mkdir $data../blob_taxonomy

blob_tax="$data../blob_taxonomy/"

for item in $blob_out
 do
  blobtools view -i $item -r all -o blob_taxonomy
  mv *blob_taxonomy* $blob_tax
  blobtools plot -i $item -r genus
  sleep 5
done

sleep 5

#filter genome assembly

mkdir $data../kept_contigs
mkdir $data../filtered_assembly
mkdir $data../final_genome

kept_con="$data../kept_contigs/"
filt_asm="$data../filtered_assembly/"
e=1
final="$data../final_genome/"

for item in $blob_tax
 do
  grep -v '#' $item | awk -F'\t' '$2 > 500'
  grep -v '#' $item | awk -F'\t' '$2 < 500'
  grep -v '#' $item | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 > 5'
  grep -v '##' $item | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 > 20' | awk -F'\t' '{print $1}' > "$i>
  assem=$(ls $spadesout*contigs.fasta | head -n"$e" | tail -n1)
  filter_contigs_by_list.py $assem contigs_to_keep_len500_cov20.txt filtered.fasta
  mv contigs_to_keep_len500_cov20.txt "$assem"_contigs_to_keep_len500_cov20.txt
  mv filtered.fasta "$assem"_filtered.fasta
  grep -f "$assem"_contigs_to_keep_len500_cov20.txt $item | awk '{w = w + $2; e = e + $5 * $2;} END>
  blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -se>
  mv "$assem"_contigs_to_keep_len500_cov20.txt $kept_con
  mv "$assem"_filtered.fasta $filt_asm
  mv *genome_vs_univec.6* $final
  e=$(($e +1))
  sleep 5
done
