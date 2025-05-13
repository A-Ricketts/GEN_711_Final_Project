#! /bin/bash

source activate genomics

data=$1
mega_out="$data../megablast_out/"
vs_contig="$data../vs_contig/"
con_db="$data../contigs_db/"
sixtS_seq="$data../16S_sequences/"
protein_ab="$data../protein_abundances/"
prokka="$data../prokka_output/"
quast="$data../quast_results/"
spadesout="$data../spades_assembly/"
trimmed_reads="$data../trimmed-reads/"
rawreads="$data../fastqc_raw-reads/"
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
univec="/home/users/amr1230/Final_Project/GEN_711_Final_Project/Univec/UniVec"

# visualize genome

mkdir $data../visual_genome

vis_gen="$data../visual_genome/"
q=1

for item in $filt_asm*
 do
  fasta=$(ls $spadesout | grep -v "contigs.fasta." | grep "contigs.fasta" | head -n"$q" | tail -n1)
  gbk=$(ls $prokka/*/*.gbk | head -n"$q" | tail -n1)
#  gbk=$(ls $prokdir/*.gbk)
  /tmp/genome_back/fix_your_gbk.py -g $gbk -f $item -o "$vis_gen$fasta".gbk
  q=$(($q + 1))
done
