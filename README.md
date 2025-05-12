# GEN_711_Final_Project

## GEN 711 Final Project

## Mystery Bacterial Genome

## Member 

Andrew Ricketts

## Background

Four mystery DNA sequences were provided that represent two different bacterial individuals, two files per bacterium. The data was provided in FASTQ files and will be processed to construct the genome of the bacteria to identify the genus and species. 

## Methods

The data was provided via the mystery bacterial genome tutorial. The files are 250bp paired-end Illumina HiSeq 2500 reads. The analysis was conducted on the RON server. 

### FastQC
FastQC runs a quality control of provided fastq files on both forward and reverse reads separately. The ouput of the program is an html file for each fastq file which can be accessed through a web browser. It provides multiple statistics and blots to show average quality and gc content of reads. 
### Trimmomatic
Trimmomatic trims off low quality bases, typically near the ends of reads, and adapter sequnces from the reads. Fastq files of both forward and reverse reads are provided as the input. The ouput is 4 fastq files: two files for the paired forward and reverse reads; and two unpaired files for both forward and reverse. 
### Spades
Spades is a genome assembler best used for bacteria. The four trimmed sequence files are inputted into the program which then aligns the reads and attempts to construct a genome. The output is multiple versions of a contigs fasta file which has the assembeled contigs. Multiple contigs are often produced when they cannot link with other contigs or there is contamination. 
### QUAST
QUAST assesses quality of the genome asembly. It looks at genome fragmentation, whether one or a few contigs is prepresentative of the genome or if there are manby small contigs. It intakes the contigs fasta file and provides quality reports in multiple file types. 
### BUSCO
BUSCO assesses the completeness of the assembled genome in the contigs fasta file. The program compares a large set of orthologous genes which 90% of all orgranisms in a certain group contain, the group in this analysis was bacteria. A complete genome would have most of these genes. The output is multiple files which report the completeness. 
### PROKKA
PROKKA is a genome annotation tool for prokaryotes which analyses the genome for where gene sequences are likely. It reports the locations and names of known genes, unknown predicted genes, and different RNAs. The contigs fasta file is inputted and multiple files are produced which each show different data like nucleotide sequence, amino acid sequence, etc.
### Extract Sequences
The extract_sequences program is utilized to find a certain gene or RNA in PROKKA ffn files. It helps extract certain sequences to compare with known organism sequences. The name of the gene or RNA is identified to find the correct sequence. The output can be put into a fasta file. 
### BLAST
#### BLAST Database
makeblastdb is a program which creates a local database of seqeunce data to quickly compare sequences with other org anisms. 
#### Local BLAST
One option for the blastn program is one against the local database. A gene or RNA seqeunce in a fasta file can be provided into the program which will compare the provided sequence with that of other organisms in the database. The output type must be provided, here a tsv file was used. Likely organism IDs are provided in the output.
#### Entire Assembly BLAST
The authors of the tutorial provided the blob_blast.sh program which BLAST's the entire genome against the local database and prepares a file which can later be used in the blob tools program. The input is the full genome fasta file and outputs a megablast.out file. 
### BWA

### Samtools

### Blobtools

### Filter Contigs

### UniVec BLAST 


## Findings

## References
