#! /bin/bash

source activate genomics

data=$1
spadesout="$data../spades_assembly/"
sixtS_seq="$data../16S_sequences/"
trimmed_reads="$data../trimmed-reads/"
con_db="$data../contigs_db/"
vs_contig="$data../vs_contig/"
mega_out="$data../megablast_out/"
raw_map="$data../raw_mapped/"
sort_map="$data../sorted_mapped/"
flag="$data../flagstat/"
bed_out="$data../bed_cov_out/"
cov_table="$data../cov_table/"
blob_out="$data../blob_out/"
blob_tax="$data../blob_taxonomy/"
kept_con="$data../kept_contigs/"
filt_asm="$data../filtered_assembly/"
final="$data../final_genome/"
e=1
univec="/home/users/amr1230/Final_Project/GEN_711_Final_Project/Univec/UniVec"

for item in $blob_tax*
 do
  grep -v '#' $item | awk -F'\t' '$2 > 500'
  grep -v '#' $item | awk -F'\t' '$2 < 500'
  grep -v '#' $item | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 > 5'
  grep -v '##' $item | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 > 20' | awk -F'\t' '{print $1}' > contigs_to_keep_len500_cov20.txt
  assem=$(ls $spadesout*contigs.fasta | head -n"$e" | tail -n1)
  mv contigs_to_keep_len500_cov20.txt "$assem"_contigs_to_keep_len500_cov20.txt
  filter_contigs_by_list.py $assem "$assem"_contigs_to_keep_len500_cov20.txt "$assem"_filtered.fasta
  grep -f "$assem"_contigs_to_keep_len500_cov20.txt $item | awk '{w = w + $2; e = e + $5 * $2;} END {print e/w}'
  blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -query "$assem"_filtered.fasta -subject $univec -outfmt 6 -out genome_vs_univec.6
  mv "$assem"_contigs_to_keep_len500_cov20.txt $kept_con
  mv "$assem"_filtered.fasta $filt_asm
  mv *genome_vs_univec.6* $final
  e=$(($e +1))
  sleep 5
done

