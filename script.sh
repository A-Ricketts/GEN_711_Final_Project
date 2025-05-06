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

rawreads="$data../fastqc_raw-reads/"

for item in $data*
 do
  fastqc $item -o $rawreads
done

sleep 5

#adpater and quality trimming

for item in $data*R1*
 do
  reverse=$(ls $data* | grep -A1 "$item" | tail -n1)
  trim_scriptV2.sh $item $data$reverse
  sleep 5
done

mv trimmed-reads "$data"../trimmed-reads

sleep 5

#genome assembly

mkdir $data../spades_assembly

spadesout="$data../spades_assembly/"
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
      mv "$spadesout"contigs.fasta "$spadesout$item"_contigs.fasta
      mv "$spadesout"spades.log "$spadesout$item"_spades.log
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
      rm -r $item
    fi
  fi
done

#genome assessment

mkdir $data../quast_results

quast="$data../quast_results/"

for item in $spadesout*contigs.fasta
 do
  quast.py $item -o $quast
  for file in $(ls $quast)
   do
    if [[ "$file" != input_* ]]
     then
      g=1
      mv $quast$file $item$file
      mv $item$file "$quast"input_$g/
      g=$(($g + 1))
    fi
  done
  sleep 5
done

sleep 5

#busco

mkdir "$data"../busco_results

buscoo="$data../busco_results/"
buscod=""$data"busco-results/"

for item in $(ls $spadesout*contigs.fasta)
 do
  for file in $spadesout*contigs.fasta
   do
    if [[ "$file" == *$item ]]
     then
      busco -i $file -m genome -o "$item"busco-results --out_path $data --force -l bacteria
      mv "$buscod"/run_bacteria_odb10/full_table.tsv "$buscod$item"_full_table.tsv
      mv "$buscod"run_bacteria_odb10/busco_sequences/single_copy_busco_sequences/ "$buscod$item"_single_copy_busco_sequences/
      for items in $buscod
       do
	if [[ "$items" != *$item* ]]
	 then
	  rm -r $items
          rm $items
	fi
      done
      cp $buscod* $buscoo
      rm -r $buscod
      rm -r "$data"busco_downloads/
    fi
  done
 sleep 5
done

sleep 5

#genome annotation

mkdir $data../prokka_output/

prokka="$data../prokka_output/"

for item in $spadesout*contigs.fasta
 do
  prokka $item --outdir $prokka --force --cpus 24 --mincontiglen 200
  sleep 5
done

sleep 5

mkdir $data../protein_abundances

protein_ab="$data../protein_abundances/"

for item in $prokka*.gff
 do
  grep -o "product=.*" $item | sed 's/product=//g' | sort | uniq -c | sort -nr > protein_abundances.txt
  mv protein_abundances.txt "$item"_protein_abundances.txt
  mv "$item"_protein_abundances.txt $protein_ab
  sleep 5
done

sleep 5

###organism identification

#16S sequence

mkdir $data../16S_sequences

sixtS_seq="$data../16S_sequences/"

for item in $prokka*.ffn
 do
  extract_sequences "16S ribosomal RNA" $item > 16S_sequence.fasta
  mv 16S_sequence.fasta "$item"_16S_sequence.fasta
  mv "$item"_16S_sequence.fasta $sixtS_seq
  sleep 5
done

sleep 5

#BLAST

mkdir $data../contigs_db

con_db="$data../contigs_db/"

for item in $(ls $spadesout*contigs.fasta)
 do
  for file in $spadesout*contigs.fasta
   do
    if [[ "$file" == *$item ]]
     then
      makeblastdb -in $file -dbtype nucl -out "$item"_contigs_db -title "$file"_title_contigs_db
      mv "$data$item"_contigs_db* $con_db
      mv "$file"_title_contigs_db $con_db
    fi
    sleep 5
done

sleep 5

mkdir $data../vs_contig

vs_contig="$data../vs_contig/"
d=1

for item in $sixtS_seq*
 do
  cont_db=$(ls $con_db* | grep "title_contigs_db" | head -n"$d" | tail -n1)
  blastn -query $item -db $cont_db -out 16S_vs_contigs_6.tsv -outfmt 6
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
  forward=$(grep "R1" $trimmed_reads | grep -v "unpaired" | head -n"$a" | tail -n1)
  reverse=$(grep "R2" $trimmed_reads | grep -v "unpaired" | head -n"$a" | tail -n1)
  bwa index $item
  bwa mem -t 24 $item $forward $reverse > "$item"_raw_mapped.sam
  mv "$item"_raw_mapped.sam $raw_map
  a=$(($a + 1))
  sleep 5
done

sleep 5

mkdir $data../sorted_mapped

sort_map="$data../sorted_mapped/"

for item in $raw_map*
 do
  samtools view -@ 24 -Sb $item | samtools sort -@ 24 -o "$item"_sorted_mapped.bam
  mv "$item"_sorted_mapped.bam $sort_map
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
  bedtools genomecov -ibam $item > "$item"_coverage.out
  mv "$item"_coverage.out $bed_out
  sleep 5
done

sleep 5

mkdir $data../cov_table

cov_table="$data../cov_table/"
b=1

for item in $bed_out*
 do
  fasta=$(head -n"$b" $spadesout*contigs.fasta | tail -n1)
  gen_input_table.py --isbedfiles $fasta $item > "$item"_coverage_table.tsv
  mv "$item"_coverage_table.tsv $cov_table
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
  sort_bam=$(head -n"$c" $sort_map* | tail -n1)
  mega=$(head -n"$c" $mega_out* | tail -n1)
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
  mv ./*blob_taxonomy* $blob_tax
  blobtools plot -i $item -r genus
  sleep 5
done
#####################################
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
  grep -v '##' $item | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 > 20' | awk -F'\t' '{print $1}' > "$item"contigs_to_keep_len500_cov20.txt
  assem=$(head -n"$e" $spadesout*contigs.fasta | tail -n1)
  filter_contigs_by_list.py $assem contigs_to_keep_len500_cov20.txt filtered.fasta
  mv contigs_to_keep_len500_cov20.txt "$item"_contigs_to_keep_len500_cov20.txt
  mv filtered.fasta "$item"_filtered.fasta
  e=$(($e +1))
  grep -f "$item"_contigs_to_keep_len500_cov20.txt $item | awk '{w = w + $2; e = e + $5 * $2;} END {print e/w}'
  blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -query "$item"_filtered.fasta -subject UniVec  -outfmt 6 -out genome_vs_univec.6
  mv "$item"_contigs_to_keep_len500_cov20.txt $kept_con
  mv "$item"_filtered.fasta $filt_asm
  mv *genome_vs_univec.6* $final
  sleep 5
done
