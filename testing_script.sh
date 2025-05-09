#! /bin/bash
source activate genomics

data=$1
spadesout="$data../spades_assembly/"
sixtS_seq="$data../16S_sequences/"
trimmed_reads="$data../trimmed-reads/"
con_db="$data../contigs_db/"
vs_contig="$data../vs_contig/"
d=1


for item in $sixtS_seq*
 do
  cont_db=$(realpath $con_db* | head -n"$d" | tail -n1)
  cona_db=$(ls $con_db | head -n"$d" | tail -n1)
  blastn -db "$cont_db"/"$cona_db"_contigs_db -query $item -outfmt 6 -out 16S_vs_contigs_6.tsv
#  while [ ! -f *16S_vs_contigs_6.tsv ]
#   do
#    echo hi
#    sleep 3
#  done
  mv 16S_vs_contigs_6.tsv "$cona_db"_16S_vs_contigs_6.tsv
  mv "$cona_db"_16S_vs_contigs_6.tsv $vs_contig
  d=$(($d + 1))
  sleep 5
done
