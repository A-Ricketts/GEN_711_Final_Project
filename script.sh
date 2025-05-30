#! /bin/bash

source activate genomics

#estimate of coverage of reads

data=$1
## data is directory with all files

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

for item in $data*fastq*
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

for item in $(ls $spadesout)
 do
  if [[ "$item" == *contigs.fasta ]]
   then
    for file in $spadesout*contigs.fasta
     do
      if [[ "$file" == *$item ]]
       then
        busco -i $file -m genome -o "$item"_busco-results --out_path $data --force -l bacteria
	buscod="$data$item""_busco-results/"
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
  fi
 sleep 5
done

sleep 5

#genome annotation

mkdir $data../prokka_output/

prokka="$data../prokka_output/"
o=1

for item in $(ls $spadesout)
 do
  if [[ "$item" == *contigs.fasta ]]
  then
    for file in $spadesout*contigs.fasta
     do
      if [[ "$file" == *$item ]]
       then
        prokka $file --outdir "$prokka$item/" --force --cpus 24 --mincontiglen 200
      fi
    done
  fi
  sleep 5
done

sleep 5

mkdir $data../protein_abundances

protein_ab="$data../protein_abundances/"
m=1

for item in $prokka*/*.gff
 do
  fasta=$(ls $spadesout*contigs.fasta | head -n"$m" | tail -n1)
  grep -o "product=.*" $item | sed 's/product=//g' | sort | uniq -c | sort -nr > protein_abundances.txt
  mv protein_abundances.txt "$fasta"_protein_abundances.txt
  mv "$fasta"_protein_abundances.txt $protein_ab
  m=$(($m + 1))
  sleep 5
done

sleep 5

###organism identification

#16S sequence

mkdir $data../16S_sequences

sixtS_seq="$data../16S_sequences/"
n=1

for item in $prokka*/*.ffn
 do
  fasta=$(ls $spadesout*contigs.fasta | head -n"$n" | tail -n1)
  extract_sequences "16S ribosomal RNA" $item > 16S_sequence.fasta
  mv 16S_sequence.fasta "$fasta"_16S_sequence.fasta
  mv "$fasta"_16S_sequence.fasta $sixtS_seq
  n=$(($n + 1))
  sleep 5
done

sleep 5

#BLAST

mkdir $data../contigs_db

con_db="$data../contigs_db/"

for item in $(ls $spadesout)
 do
  if [[ "$item" == *contigs.fasta ]]
  then
    for file in $spadesout*contigs.fasta
     do
      if [[ "$file" == *$item ]]
       then
	mkdir $con_db$item
        makeblastdb -in $file -dbtype nucl -out "$con_db$item"/"$item"_contigs_db
      fi
    done
  fi
  sleep 5
done

sleep 5

mkdir $data../vs_contig

vs_contig="$data../vs_contig/"
d=1

for item in $sixtS_seq*
 do
  cont_db=$(realpath $con_db* | head -n"$d" | tail -n1)
  cona_db=$(ls $con_db | head -n"$d" | tail -n1)
  blastn -db "$cont_db"/"$cona_db"_contigs_db -query "$item" -outfmt 6 -out 16S_vs_contigs_6.tsv
  mv 16S_vs_contigs_6.tsv "$cona_db"_16S_vs_contigs_6.tsv
  mv "$cona_db"_16S_vs_contigs_6.tsv $vs_contig
  d=$(($d + 1))
  sleep 5
done

sleep 5

mkdir $data../megablast_out

mega_out="$data../megablast_out/"

for item in $spadesout*contigs.fasta
 do
  blast-ncbi-nt.sh $item
  while [ ! -f *megablast.out ]
   do
    sleep 30
  done
  mv *megablast.out $mega_out
done

sleep 5

#read mapping

mkdir $data../raw_mapped
raw_map="$data../raw_mapped/"
a=1

for item in $spadesout*contigs.fasta
 do
  forward=$(ls $trimmed_reads* | grep "R1" | grep -v "unpaired" | head -n"$a" | tail -n1)
  reverse=$(ls $trimmed_reads* | grep "R2" | grep -v "unpaired" | head -n"$a" | tail -n1)
  bwa index $item
  bwa mem -t 24 $item $forward $reverse > raw_mapped.sam
  fasta=$(ls $spadesout*contigs.fasta | head -n"$a" | tail -n1)
  mv raw_mapped.sam "$fasta"_raw_mapped.sam
  mv "$fasta"_raw_mapped.sam $raw_map
  a=$(($a + 1))
  sleep 5
done

sleep 5

mkdir $data../sorted_mapped

sort_map="$data../sorted_mapped/"
i=1

for item in $raw_map*
 do
  samtools view -@ 24 -Sb $item | samtools sort -@ 24 -o sorted_mapped.bam
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
h=1

for item in $sort_map*bam
 do
  samtools flagstat $item > "$item"_flagstat.txt
  mv "$item"_flagstat.txt $flag
  samtools index $item
  bedtools genomecov -ibam $item > coverage.out
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
  /tmp/genome_back/gen_input_table.py --isbedfiles $fasta $item > coverage_table.tsv
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
  sort_bam=$(ls $sort_map*.bam | head -n"$c" | tail -n1)
  mega=$(ls $mega_out* | head -n"$c" | tail -n1)
  fasta=$(ls $spadesout | grep -v "contigs.fasta." | grep "contigs.fasta" | head -n"$c" | tail -n1)
  blobtools create -i $item -b $sort_bam -t $mega -o "$blob_out$fasta"_blob_out
  c=$(($c + 1))
  sleep 5
done

sleep 5

mkdir $data../blob_taxonomy

blob_tax="$data../blob_taxonomy/"
m=1

for item in $blob_out*json
 do
  fasta=$(ls $spadesout | grep -v "contigs.fasta." | grep "contigs.fasta" | head -n"$m" | tail -n1)
  blobtools view -i $item -r all -o "$blob_tax$fasta"_blob_taxonomy
  blobtools plot -i $item -r genus
  mv ./*blobplot* $blob_tax
  m=$(($m + 1))
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
univec="/home/users/amr1230/Final_Project/GEN_711_Final_Project/Univec/UniVec"

for item in $blob_tax*table.txt
 do
  grep -v '##' $item | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 > 20' | awk -F'\t' '{print $1}' > contigs_to_keep_len500_cov20.txt
  assem=$(ls $spadesout*contigs.fasta | head -n"$e" | tail -n1)
  mv contigs_to_keep_len500_cov20.txt "$assem"_contigs_to_keep_len500_cov20.txt
  filter_contigs_by_list.py $assem "$assem"_contigs_to_keep_len500_cov20.txt "$assem"_filtered.fasta
  grep -f "$assem"_contigs_to_keep_len500_cov20.txt $item | awk '{w = w + $2; e = e + $5 * $2;} END {print e/w}'
  blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -query "$assem"_filtered.fasta -subject $univec -outfmt 6 -out genome_vs_univec.6
  mv "$assem"_contigs_to_keep_len500_cov20.txt $kept_con
  mv "$assem"_filtered.fasta $filt_asm
  mv genome_vs_univec.6 "$assem"_genome_vs_univec.6
  mv "$assem"_genome_vs_univec.6 $final
  e=$(($e +1))
  sleep 5
done

# visualize genome

mkdir $data../visual_genome

vis_gen="$data../visual_genome/"
q=1

for item in $filt_asm*
 do
  fasta=$(ls $spadesout | grep -v "contigs.fasta." | grep "contigs.fasta" | head -n"$q" | tail -n1)
  gbk=$(ls $prokka/*/*.gbk | head -n"$q" | tail -n1)
  /tmp/genome_back/fix_your_gbk.py -g $gbk -f $item -o "$vis_gen$fasta".gbk
  q=$(($q + 1))
done
